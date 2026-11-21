"""Exporta resultados disponibles del cierre, con probabilidad exacta o cotas."""
from __future__ import annotations

import csv
from fractions import Fraction
import gzip
import json
import math

import matplotlib.pyplot as plt

from cierre_plan1 import ROOT, SALIDA, TABLAS, SPEC, bateria


def guardar_csv(nombre, filas, columnas=None):
    if not filas:
        return
    with (TABLAS / nombre).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columnas or list(filas[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(filas)


def fila_resumen(variante, filas, matriz):
    pesos = {e: Fraction(v[3]).limit_denominator(248) for e, v in matriz.items()}
    assert sum(pesos.values()) == 1
    pf = sum((pesos[f["escenario"]] for f in filas if f["falla_conjunta"]), Fraction())
    no_evaluados = set(matriz) - {f["escenario"] for f in filas}
    superior = pf + sum(pesos[e] for e in no_evaluados)
    pico_t = max(filas, key=lambda f: f["M_T_c"])
    pico_delta = max(filas, key=lambda f: f["M_delta_c"])
    return dict(variante=variante, escenarios=len(filas), D_h=pico_t["D_h"],
                max_T_336_c=max(f["MT_336_c"] for f in filas),
                max_delta_336_c=max(f["Mdelta_336_c"] for f in filas),
                max_T_c=pico_t["M_T_c"], escenario_T=pico_t["escenario"], t_T_h=pico_t["t_M_T_h"],
                max_delta_c=pico_delta["M_delta_c"], escenario_delta=pico_delta["escenario"],
                t_delta_h=pico_delta["t_M_delta_h"],
                exceden_T=sum(f["falla_T"] for f in filas),
                exceden_delta=sum(f["falla_delta"] for f in filas),
                exceden_union=sum(f["falla_conjunta"] for f in filas),
                Pf=str(pf) if not no_evaluados else "",
                Pf_cota_inferior=str(pf), Pf_cota_superior=str(superior),
                clasificacion=("excede_5pct" if pf > Fraction(1, 20) else
                               "cumple_5pct" if superior <= Fraction(1, 20) else "indeterminada"),
                max_balance_relativo=max(f["balance_relativo"] for f in filas),
                max_dt_sobre_estable=max(f["dt_sobre_estable_max"] for f in filas),
                costo_total_s=sum(f["segundos"] for f in filas))


def main():
    TABLAS.mkdir(parents=True, exist_ok=True)
    datos = []
    for ruta in sorted(SALIDA.glob("*/*.json.gz")):
        with gzip.open(ruta, "rt") as f:
            datos.append(json.load(f))
    matriz = bateria._contexto()["matriz"]
    filas = [d["fila"] for d in datos]
    guardar_csv("escenarios.csv", filas)
    resumen = [fila_resumen(v, grupo, matriz) for v in SPEC["variantes"]
               if (grupo := [f for f in filas if f["variante"] == v])]
    guardar_csv("resumen.csv", resumen)
    ejecucion = [dict(variante=d["fila"]["variante"], escenario=d["fila"]["escenario"],
                      firma=d["firma"], rss_max_kib=d["rss_max_kib"], python=d["python"],
                      numpy=d["numpy"], plataforma=d["plataforma"],
                      retiros=json.dumps(d["retiros"], ensure_ascii=False)) for d in datos]
    guardar_csv("ejecuciones.csv", ejecucion)
    ref = {f["escenario"]: f for f in filas if f["variante"] == "exacto"}
    historico = ROOT / "reports/Defensa/bateria7/tablas/resultados_por_escenario.csv"
    with historico.open(encoding="utf-8", newline="") as f:
        original = {r["escenario"]: r for r in csv.DictReader(f) if r["caso"].startswith("4 ")}
    control = []
    for e, f in ref.items():
        previo = original[e]
        dt = f["MT_336_c"] - float(previo["M_T_c"])
        dd = f["Mdelta_336_c"] - float(previo["M_delta_c"])
        control.append(dict(escenario=e, publicado_MT_c=previo["M_T_c"],
                            reproducido_MT_c=f["MT_336_c"], diferencia_MT_c=dt,
                            publicado_delta_c=previo["M_delta_c"],
                            reproducido_delta_c=f["Mdelta_336_c"], diferencia_delta_c=dd,
                            reproduce_a_3_decimales=max(abs(dt), abs(dd)) <= 0.0005))
    guardar_csv("control_336h.csv", control)
    sensibilidad = []
    for f in filas:
        if f["variante"] in ("espacio", "tiempo") and f["escenario"] in ref:
            base = ref[f["escenario"]]
            sensibilidad.append(dict(variante=f["variante"], escenario=f["escenario"],
                                     delta_MT_336_c=f["MT_336_c"] - base["MT_336_c"],
                                     delta_Mdelta_336_c=f["Mdelta_336_c"] - base["Mdelta_336_c"],
                                     delta_MT_672_c=f["M_T_c"] - base["M_T_c"],
                                     delta_Mdelta_672_c=f["M_delta_c"] - base["M_delta_c"],
                                     cambia_falla=f["falla_conjunta"] != base["falla_conjunta"]))
    guardar_csv("sensibilidad.csv", sensibilidad)
    series = []
    diagnostico = []
    fig, ejes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
    for d in datos:
        f = d["fila"]
        critico = f["escenario"] == "M4_SS__GEFS_20261116_frio"
        s = d["serie"]
        retiro_h = next(r[0] / 3600 for r in d["retiros"] if any("manta_PUR" in n for n in r[1]))
        alrededor = [r["diferencia_maxima_S25_c"] for r in s if r["tiempo_h"] == retiro_h]
        ultimas = [r["diferencia_maxima_S25_c"] for r in s
                   if r["tiempo_h"] > f["horizonte_h"] - 24]
        diagnostico.append(dict(variante=f["variante"], escenario=f["escenario"],
                                primer_retiro_h=retiro_h, diferencial_antes_c=alrededor[0],
                                diferencial_despues_c=alrededor[-1], diferencial_final_c=ultimas[-1],
                                max_delta_ultimas_24h_c=max(ultimas),
                                pico_delta_ultimas_24h=f["t_M_delta_h"] > f["horizonte_h"] - 24,
                                pico_T_ultimas_24h=f["t_M_T_h"] > f["horizonte_h"] - 24))
        for r in s:
            t = r["tiempo_h"]
            if math.isclose(t, round(t), abs_tol=1e-8) or t in (f["t_M_T_h"], f["t_M_delta_h"]):
                series.append(dict(variante=f["variante"], escenario=f["escenario"],
                                   h=t, Tmax_c=round(r["temperatura_maxima_c"], 6),
                                   S25min_c=round(r["temperatura_minima_S25_c"], 6),
                                   delta_c=round(r["diferencia_maxima_S25_c"], 6)))
        if critico:
            hs = [r["tiempo_h"] for r in s]
            for eje, col in zip(ejes, ("temperatura_maxima_c", "diferencia_maxima_S25_c")):
                eje.plot(hs, [r[col] for r in s], label=f["variante"], linewidth=1)
    guardar_csv("diagnostico_retiro.csv", diagnostico)
    # Los eventos pueden tener estados simultáneos antes/después.
    guardar_csv("curvas.csv", list({tuple(r.values()): r for r in series}.values()))
    for eje, umbral, etiqueta in zip(ejes, (70, 20), ("T máxima [°C]", "Diferencial S25 [°C]")):
        eje.axhline(umbral, color="firebrick", linestyle=":", label="Umbral")
        eje.axvline(336, color="0.4", linestyle="--", label="Cierre original")
        eje.axvline(375, color="0.5", linestyle=":", label="Retiro original")
        eje.set_ylabel(etiqueta)
        eje.grid(alpha=0.2)
    ejes[0].legend(fontsize="small", ncols=3)
    ejes[0].set_title("M4_SS · GEFS frío")
    ejes[-1].set_xlabel("Horas desde inicio de obra")
    fig.savefig(TABLAS / "curvas_criticas.png", dpi=130)


if __name__ == "__main__":
    main()
