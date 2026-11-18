"""Consolida la batería de 7 casos: tablas, comprobación de reutilización y figuras.

    python scripts/bateria7_analisis.py

Lee results/bateria7/cache y las matrices ya publicadas (reports/Semana Final/evidencia_semana_final/
comparacion_100m3h/resultados/resultados_por_escenario.csv) y escribe results/bateria7/tablas|figuras.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
PUB = RAIZ / "reports/Semana Final/evidencia_semana_final/comparacion_100m3h/resultados/resultados_por_escenario.csv"
CASOS = [("bat1_base_original", "1 Base original", "caso_ifc_original_sin_manta"),
         ("bat2_base_simultaneo", "2 Base simultáneo", None),
         ("bat3_base_largo_24h", "3 Base largo 24 h", None),
         ("bat4_plan1_original_10h", "4 Plan 1 original 10 h", "propuesto_10h_60mm_T10"),
         ("bat5_plan_simultaneo", "5 Plan simultáneo", "vertida_unica_100m3h"),
         ("bat6_plan_largo_24h", "6 Plan largo 24 h", None),
         ("bat7_control_sin_manta", "7 Control sin manta", "control_10h_T10_sin_manta")]
CURVAS = ["M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"]
CLIMAS = {"frio": 11 / 31, "central": 10 / 31, "calido": 10 / 31}
ESC = lambda c, w: f"{c}__GEFS_20261116_{w}"


def cargar(plan_id):
    d = {}
    for ruta in (SAL / "cache").glob(f"{plan_id}__*.json"):
        r = json.loads(ruta.read_text(encoding="utf-8"))
        d[r["escenario_id"]] = r
    return d


def publicados():
    out = {}
    with PUB.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[(r["plan_id"], r["escenario_id"])] = r
    return out


def main():
    (SAL / "tablas").mkdir(exist_ok=True)
    pub = publicados()
    # 1) comprobación de reutilización
    filas = []
    for pid, nombre, orig in CASOS:
        if orig is None:
            continue
        for esc, r in cargar(pid).items():
            p = pub[(orig, esc)]
            filas.append({"caso": nombre, "escenario": esc,
                          "dMT_c": r["metricas"]["M_T_c"] - float(p["M_T_c"]),
                          "dMdelta_c": r["metricas"]["M_delta_c"] - float(p["M_delta_c"]),
                          "dt_MT_h": r["metricas"]["t_M_T_h"] - float(p["t_M_T_h"]),
                          "dt_Mdelta_h": r["metricas"]["t_M_delta_h"] - float(p["t_M_delta_h"])})
    with (SAL / "tablas/comprobacion_reutilizacion.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)
    print("max |dMT|", max(abs(x["dMT_c"]) for x in filas), " max |dMdelta|", max(abs(x["dMdelta_c"]) for x in filas), f"({len(filas)} comparaciones)")

    # 2) tabla por escenario (24 por caso): casos nuevos desde cache; reutilizados desde publicado
    todo = {}
    for pid, nombre, orig in CASOS:
        for c in CURVAS:
            for w in CLIMAS:
                e = ESC(c, w)
                if orig:
                    p = pub[(orig, e)]
                    todo[(pid, e)] = dict(M_T=float(p["M_T_c"]), M_d=float(p["M_delta_c"]), t_T=float(p["t_M_T_h"]),
                                          t_d=float(p["t_M_delta_h"]), xT=(float(p["x_M_T_m"]), float(p["y_M_T_m"]), float(p["z_M_T_m"])),
                                          xS=(float(p["x_S25_min_m"]), float(p["y_S25_min_m"]), float(p["z_S25_min_m"])), fuente="publicado")
                else:
                    r = cargar(pid)[e]["metricas"]
                    todo[(pid, e)] = dict(M_T=r["M_T_c"], M_d=r["M_delta_c"], t_T=r["t_M_T_h"], t_d=r["t_M_delta_h"],
                                          xT=tuple(r["xyz_M_T_m"]), xS=tuple(r["xyz_S25_min_m"]), fuente="nuevo")
    with (SAL / "tablas/resultados_por_escenario.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["caso", "escenario", "peso", "M_T_c", "M_delta_c", "t_M_T_h", "t_M_delta_h", "x_MT", "y_MT", "z_MT", "x_S25", "y_S25", "z_S25", "falla_T", "falla_delta", "fuente"])
        for pid, nombre, _ in CASOS:
            for c in CURVAS:
                for wk, wp in CLIMAS.items():
                    e = ESC(c, wk); d = todo[(pid, e)]
                    w.writerow([nombre, e, wp / 8, f"{d['M_T']:.3f}", f"{d['M_d']:.3f}", f"{d['t_T']:.1f}", f"{d['t_d']:.1f}",
                                *[f"{v:.3f}" for v in d["xT"]], *[f"{v:.3f}" for v in d["xS"]], int(d["M_T"] > 70), int(d["M_d"] > 20), d["fuente"]])

    # 3) resumen por caso con 4 interpretaciones de P_f
    res = []
    for pid, nombre, _ in CASOS:
        def P(escs):
            tot = sum(CLIMAS[e.split("_")[-1]] for e in escs)
            return sum(CLIMAS[e.split("_")[-1]] for e in escs if todo[(pid, e)]["M_T"] > 70 or todo[(pid, e)]["M_d"] > 20) / tot
        todos = [ESC(c, w) for c in CURVAS for w in CLIMAS]
        ss4 = [ESC(c, w) for c in CURVAS[:4] for w in CLIMAS]
        nom = [ESC("M1_SS", w) for w in CLIMAS]
        peorT = max(todos, key=lambda e: todo[(pid, e)]["M_T"]); peorD = max(todos, key=lambda e: todo[(pid, e)]["M_d"])
        nF = sum(1 for e in todos if todo[(pid, e)]["M_T"] > 70 or todo[(pid, e)]["M_d"] > 20)
        res.append({"caso": nombre, "Pf_8curvas_3climas": round(P(todos), 4), "Pf_4SS": round(P(ss4), 4), "Pf_nominal_M1_SS": round(P(nom), 4),
                    "escenarios_fallidos_de_24": nF,
                    "MT_max_c": round(todo[(pid, peorT)]["M_T"], 3), "escenario_MT": peorT, "t_MT_h": todo[(pid, peorT)]["t_T"], "xyz_MT": tuple(round(v, 2) for v in todo[(pid, peorT)]["xT"]),
                    "Mdelta_max_c": round(todo[(pid, peorD)]["M_d"], 3), "escenario_Mdelta": peorD, "t_Mdelta_h": todo[(pid, peorD)]["t_d"], "xyz_S25": tuple(round(v, 2) for v in todo[(pid, peorD)]["xS"]),
                    "Mdelta_min_c": round(min(todo[(pid, e)]["M_d"] for e in todos), 3),
                    "pico_en_ultimas_24h": int(todo[(pid, peorD)]["t_d"] >= 336 - 24)})
    with (SAL / "tablas/resumen_por_caso.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    for r in res:
        print(r)


if __name__ == "__main__":
    main()
