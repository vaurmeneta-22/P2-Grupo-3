"""Tablas del candidato con aislación inferior de 100 mm (56 d y retiro hasta 150 d)."""
import csv, glob, json, sys
from pathlib import Path
import numpy as np
RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts")); sys.path.insert(0, str(RAIZ / "src"))
import defensa_busqueda as db
from defensa_consolidar import INTERP
from mcplan.thermal.decision import ResultadoEscenario, agregar_probabilidades
T = RAIZ / "reports/Defensa/tablas"
RET = [1344., 1680., 2016., 2352.]
def cargar(prefijo):
    return [json.load(open(f)) for f in glob.glob(str(RAIZ / f"results/defensa/cache2/{prefijo}*.json"))]
def filas(prefijo, retiro):
    out = []
    for d in cargar(prefijo):
        suelo = "0,8 m" if "L0_8" in d["nombre"] else "2 m"
        h = np.array(d["serie"]["h"]); dl = np.array(d["serie"]["delta_c"]); tm = np.array(d["serie"]["Tmax_c"])
        r = dict(suelo=suelo, curva=d["curva"], clima=d["clima"], peso_clima=db.CLIMAS[d["clima"]][1], D_h=d["D_h"])
        for k, e in (("14d", d["h14"]), ("fin", d["hfin"])):
            r.update({f"M_T_{k}": e["M_T_c"], f"t_MT_{k}_h": e["t_M_T_h"], f"xyz_MT_{k}": e["xyz_M_T_m"],
                      f"M_d_{k}": e["M_delta_c"], f"t_Md_{k}_h": e["t_M_delta_h"], f"xyz_S25_{k}": e["xyz_S25_min_m"]})
        m = h <= 672; r["M_T_28d"] = tm[m].max(); r["M_d_28d"] = dl[m].max()
        m = h <= 1344; r["M_T_56d"] = tm[m].max(); r["M_d_56d"] = dl[m].max(); r["t_Md_56d_h"] = h[m][dl[m].argmax()]
        if retiro:
            b = RET + [3600.1]
            for j in range(4):
                mm = (h >= b[j]) & (h < b[j + 1]); r[f"Md_post_retiro{j+1}"] = dl[mm].max()
            mm = h >= 2352; r["Md_post_retiro_completo"] = dl[mm].max(); r["t_post_completo_h"] = h[mm][dl[mm].argmax()]
            last = h >= 3600 - 240
            r["pend_dT_10d"] = np.polyfit(h[last] / 24, dl[last], 1)[0]; r["Tmax_final"] = tm[-1]; r["delta_final"] = dl[-1]
        out.append(r)
    return sorted(out, key=lambda r: (r["suelo"], r["curva"], r["clima"]))
def planes(rows, horizontes):
    res = []
    for suelo in ("2 m", "0,8 m"):
        R = [r for r in rows if r["suelo"] == suelo]
        for ip, pc in INTERP.items():
            fila = dict(suelo=suelo, interpretacion=ip, D_h=R[0]["D_h"])
            for hz in horizontes:
                lst = [ResultadoEscenario("p", f"{r['curva']}__{r['clima']}", pc[r["curva"]] * r["peso_clima"], r[f"M_T_{hz}"], r[f"M_d_{hz}"]) for r in R if r["curva"] in pc]
                a = agregar_probabilidades(lst)[0]
                fila.update({f"P_T_{hz}": a["P_T"], f"P_delta_{hz}": a["P_delta"], f"P_f_{hz}": a["P_f"], f"maxMd_{hz}": max(x.M_delta_c for x in lst), f"maxMT_{hz}": max(x.M_T_c for x in lst), f"margen_{hz}": 20 - max(x.M_delta_c for x in lst)})
            res.append(fila)
    return res
def guardar(nombre, rows):
    with open(T / nombre, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
a = filas("ai100_", False); b = filas("ai100ret_", True)
guardar("ai100_escenarios_56d.csv", a); guardar("ai100_retiro_escenarios_150d.csv", b)
pa = planes(a, ("14d", "28d", "56d")); pb = planes(b, ("56d", "fin"))
guardar("ai100_planes_56d.csv", pa); guardar("ai100_retiro_planes_150d.csv", pb)
for r in pa + pb: print({k: (round(float(v), 3) if isinstance(v, (float,)) else v) for k, v in r.items()})
for nombre, rows, k in (("56d", a, "56d"), ("150d", b, "fin")):
    for suelo in ("2 m", "0,8 m"):
        R = [r for r in rows if r["suelo"] == suelo]; p = max(R, key=lambda r: r[f"M_d_{k}"])
        print(nombre, suelo, "peor", p["curva"], p["clima"], round(p[f"M_d_{k}"], 3), round(p[f"t_Md_{k}_h"], 1), [round(x, 3) for x in p.get(f"xyz_S25_{k}", [])],
              "Tmax", round(max(r[f"M_T_{k}"] for r in R), 3), "| top5:", [(r["curva"], r["clima"], round(r[f"M_d_{k}"], 2)) for r in sorted(R, key=lambda r: -r[f"M_d_{k}"])[:3]])
for suelo in ("2 m", "0,8 m"):
    R = [r for r in b if r["suelo"] == suelo]
    print("retiro", suelo, "max post r1..4", [round(max(r[f"Md_post_retiro{j}"] for r in R), 2) for j in (1, 2, 3, 4)], "completo", round(max(r["Md_post_retiro_completo"] for r in R), 3),
          "pend", round(max(r["pend_dT_10d"] for r in R), 3), round(min(r["pend_dT_10d"] for r in R), 3), "delta_final max", round(max(r["delta_final"] for r in R), 2), "Tmax_final", round(min(r["Tmax_final"] for r in R), 1), round(max(r["Tmax_final"] for r in R), 1))
    print(" exceden 56d:", sum(1 for r in R if r["M_d_fin"] > r["M_d_56d"] + 1e-6), "max exceso", round(max(r["M_d_fin"] - r["M_d_56d"] for r in R), 2))
