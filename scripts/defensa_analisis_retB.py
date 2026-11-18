"""Matriz completa de la regla B (24 escenarios × 2 suelos, 3400 h): tablas y probabilidades."""
import csv, glob, json, sys
from pathlib import Path
import numpy as np
RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts")); sys.path.insert(0, str(RAIZ / "src"))
import defensa_busqueda as db
from defensa_consolidar import INTERP
from mcplan.thermal.decision import ResultadoEscenario, agregar_probabilidades
T = RAIZ / "reports/Defensa/tablas"
RET = [24.0 * d for d in (56, 77, 98, 119)]
rows = []
for f in glob.glob(str(RAIZ / "results/defensa/cache2/retB_*.json")):
    d = json.load(open(f)); suelo = "0,8 m" if "L0_8" in d["nombre"] else "2 m"
    h = np.array(d["serie"]["h"]); dl = np.array(d["serie"]["delta_c"]); tm = np.array(d["serie"]["Tmax_c"])
    r = dict(suelo=suelo, curva=d["curva"], clima=d["clima"], peso_clima=db.CLIMAS[d["clima"]][1], D_h=d["D_h"])
    for k, e in (("14d", d["h14"]), ("fin", d["hfin"])):
        r.update({f"M_T_{k}": e["M_T_c"], f"t_MT_{k}_h": e["t_M_T_h"], f"xyz_MT_{k}": e["xyz_M_T_m"],
                  f"M_d_{k}": e["M_delta_c"], f"t_Md_{k}_h": e["t_M_delta_h"], f"xyz_S25_{k}": e["xyz_S25_min_m"]})
    for etq, lim in (("28d", 672), ("56d", 1344)):
        m = h <= lim; r[f"M_T_{etq}"] = tm[m].max(); r[f"M_d_{etq}"] = dl[m].max()
    r["margen_20_C"] = 20 - r["M_d_fin"]
    b = RET + [1e9]
    for j in range(4):
        m = (h >= b[j]) & (h < b[j + 1]); r[f"Md_post_retiro{j+1}"] = dl[m].max()
    r["Md_pre_retiro"] = dl[h < RET[0]].max()
    last = h >= h[-1] - 240; r["pend_dT_10d_por_dia"] = np.polyfit(h[last] / 24, dl[last], 1)[0]
    r["delta_final"] = dl[-1]; r["Tmax_final"] = tm[-1]; r["t_fin_d"] = h[-1] / 24
    r["pico_cerca_del_final"] = int(r["t_Md_fin_h"] > h[-1] - 240)
    rows.append(r)
rows.sort(key=lambda r: (r["suelo"], r["curva"], r["clima"]))
def guardar(n, R):
    with open(T / n, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(R[0])); w.writeheader(); w.writerows(R)
guardar("retB_escenarios.csv", rows)
out = []
for suelo in ("2 m", "0,8 m"):
    R = [r for r in rows if r["suelo"] == suelo]
    for ip, pc in INTERP.items():
        fila = dict(suelo=suelo, interpretacion=ip, D_h=R[0]["D_h"])
        for hz in ("14d", "28d", "56d", "fin"):
            L = [ResultadoEscenario("p", f"{r['curva']}__{r['clima']}", pc[r["curva"]] * r["peso_clima"], r[f"M_T_{hz}"], r[f"M_d_{hz}"]) for r in R if r["curva"] in pc]
            a = agregar_probabilidades(L)[0]
            fila.update({f"P_T_{hz}": a["P_T"], f"P_delta_{hz}": a["P_delta"], f"P_f_{hz}": a["P_f"], f"maxMd_{hz}": max(x.M_delta_c for x in L), f"maxMT_{hz}": max(x.M_T_c for x in L)})
        fila["margen_fin"] = 20 - fila["maxMd_fin"]
        out.append(fila); print({k: (round(float(v), 3) if isinstance(v, float) else v) for k, v in fila.items()})
guardar("retB_planes.csv", out)
for suelo in ("2 m", "0,8 m"):
    R = [r for r in rows if r["suelo"] == suelo]
    top = sorted(R, key=lambda r: -r["M_d_fin"])[:4]
    print(suelo, "top:", [(r["curva"], r["clima"], round(r["M_d_fin"], 2), round(r["t_Md_fin_h"] / 24, 1), [round(x, 2) for x in r["xyz_S25_fin"]]) for r in top])
    print("  fallos:", sum(1 for r in R if r["M_d_fin"] > 20 or r["M_T_fin"] > 70), "| peor margen", round(min(r["margen_20_C"] for r in R), 2), "| Tmax max", round(max(r["M_T_fin"] for r in R), 2))
    print("  picos tras retiro k (max sobre escenarios):", [round(max(r[f"Md_post_retiro{j}"] for r in R), 2) for j in (1, 2, 3, 4)], "pre", round(max(r["Md_pre_retiro"] for r in R), 2))
    print("  pend final max/min", round(max(r["pend_dT_10d_por_dia"] for r in R), 3), round(min(r["pend_dT_10d_por_dia"] for r in R), 3), "| delta_final max", round(max(r["delta_final"] for r in R), 2), "| picos cerca del final:", sum(r["pico_cerca_del_final"] for r in R), "| instantes pico dias min/max", round(min(r["t_Md_fin_h"] for r in R) / 24, 1), round(max(r["t_Md_fin_h"] for r in R) / 24, 1))
    print("  max controlado por retiro1 en", sum(1 for r in R if abs(r["t_Md_fin_h"] / 24 - 56) < 8), "de", len(R), "; por base pre-retiro", sum(1 for r in R if r["t_Md_fin_h"] < RET[0]))
