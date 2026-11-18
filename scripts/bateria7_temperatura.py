"""Sensibilidad a la temperatura de colocación del plan 1 (10/11/12/13 °C)."""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
PUB = RAIZ / "reports/Semana Final/evidencia_semana_final/comparacion_100m3h/resultados/resultados_por_escenario.csv"
CASOS = {10: ("propuesto_10h_60mm_T10", "bat4_plan1_original_10h"), 11: (None, "bat4_T11"), 12: (None, "bat4_T12"), 13: (None, "bat4_T13")}
CURVAS = ["M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"]
WC = {"frio": 0.354838709677, "central": 0.322580645161, "calido": 0.322580645161}
pub = {}
for r in csv.DictReader(PUB.open(encoding="utf-8")):
    pub[(r["plan_id"], r["escenario_id"])] = r
def esc(c, w): return f"{c}__GEFS_20261116_{w}"
res = {}
for T, (orig, bat) in CASOS.items():
    d = {}
    for c in CURVAS:
        for w in WC:
            e = esc(c, w)
            if orig:
                r = pub[(orig, e)]; d[e] = (float(r["M_T_c"]), float(r["M_delta_c"]), float(r["t_M_T_h"]), float(r["t_M_delta_h"]))
            else:
                m = json.loads((SAL / "cache" / f"{bat}__{e}.json").read_text(encoding="utf-8"))["metricas"]
                d[e] = (m["M_T_c"], m["M_delta_c"], m["t_M_T_h"], m["t_M_delta_h"])
    res[T] = d
# comprobación: los 3 escenarios del caso 10 °C recalculado coinciden con el publicado
for e in ("M1_SS__GEFS_20261116_central", "M4_SS__GEFS_20261116_frio", "M4_SS__GEFS_20261116_central"):
    m = json.loads((SAL / "cache" / f"bat4_plan1_original_10h__{e}.json").read_text(encoding="utf-8"))["metricas"]
    assert abs(m["M_delta_c"] - res[10][e][1]) < 1e-9 and abs(m["M_T_c"] - res[10][e][0]) < 1e-9
filas = []
for T, d in res.items():
    pT = sum(WC[e.split("_")[-1]] / 8 for e in d if d[e][0] > 70)
    pD = sum(WC[e.split("_")[-1]] / 8 for e in d if d[e][1] > 20)
    pF = sum(WC[e.split("_")[-1]] / 8 for e in d if d[e][0] > 70 or d[e][1] > 20)
    eD = max(d, key=lambda e: d[e][1]); eT = max(d, key=lambda e: d[e][0])
    nf = sum(1 for e in d if d[e][1] > 20)
    ser = json.loads((SAL / "cache" / f"{CASOS[T][1]}__{eD}.json").read_text(encoding="utf-8"))["serie"] if (SAL / "cache" / f"{CASOS[T][1]}__{eD}.json").exists() else None
    pend = None
    if ser:
        t = np.asarray(ser["tiempo_h"]); y = np.asarray(ser["diferencia_maxima_S25_c"])
        pend = (y[-1] - y[t >= t[-1] - 24][0]) / 24.0
    filas.append(dict(T_col_c=T, M_T_max_c=round(d[eT][0], 3), esc_M_T=eT, t_M_T_h=round(d[eT][2], 1), M_delta_max_c=round(d[eD][1], 3),
                      margen_delta_c=round(20 - d[eD][1], 3), margen_T_c=round(70 - d[eT][0], 3), esc_M_delta=eD, t_M_delta_h=round(d[eD][3], 1),
                      escenarios_fallan_delta=nf, P_T=pT, P_delta=pD, P_f=pF, pendiente_final_c_por_h=None if pend is None else round(pend, 4)))
with (SAL / "tablas/sensibilidad_temperatura_plan1.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(filas[0])); w.writeheader(); w.writerows(filas)
for r in filas: print(r)
print("por grado:")
for a, b in zip(filas, filas[1:]):
    print(f"{a['T_col_c']}->{b['T_col_c']}: dTmax {b['M_T_max_c']-a['M_T_max_c']:+.3f}  dDelta {b['M_delta_max_c']-a['M_delta_max_c']:+.3f}  crit {a['esc_M_delta']==b['esc_M_delta']}")
# por escenario: ΔM_Δ por °C
dd = {e: (res[13][e][1] - res[10][e][1]) / 3 for e in res[10]}; dT = {e: (res[13][e][0] - res[10][e][0]) / 3 for e in res[10]}
print("dMdelta/°C por escenario: min %.3f max %.3f; dMT/°C min %.3f max %.3f" % (min(dd.values()), max(dd.values()), min(dT.values()), max(dT.values())))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
Ts = [r["T_col_c"] for r in filas]
ax[0].plot(Ts, [r["M_T_max_c"] for r in filas], "o-", color="#d62728"); ax[0].axhline(70, color="k"); ax[0].text(10.05, 70.6, "límite 70 °C", fontsize=9); ax[0].set_title("Tmax máxima de los 24 escenarios"); ax[0].set_ylim(40, 75)
ax[1].plot(Ts, [r["M_delta_max_c"] for r in filas], "o-", color="#1f77b4"); ax[1].axhline(20, color="k"); ax[1].text(10.05, 20.03, "límite 20 °C", fontsize=9); ax[1].set_title("ΔT máxima de los 24 escenarios")
ax[2].plot(Ts, [100 * r["P_f"] for r in filas], "o-", color="#2ca02c"); ax[2].axhline(5, color="k"); ax[2].text(10.05, 5.25, "límite 5 %", fontsize=9); ax[2].set_title("Pf conjunta")
for a in ax: a.set_xlabel("Temperatura de colocación uniforme [°C]"); a.set_xticks(Ts); a.grid(alpha=.3)
ax[0].set_ylabel("Tmax máxima de los 24 escenarios [°C]"); ax[1].set_ylabel("Diferencia térmica ΔT máxima [°C]"); ax[2].set_ylabel("Probabilidad de fallo Pf [%]")
fig.suptitle("Plan 1 (10 h, PUR 60 mm desde 15 h): sensibilidad a la temperatura de colocación; 336 h"); fig.tight_layout()
fig.savefig(SAL / "figuras/sensibilidad_T_colocacion.png", dpi=130); plt.close(fig)
cols = {10: "#1f77b4", 11: "#2ca02c", 12: "#ff7f0e", 13: "#d62728"}
for e in ("M4_SS__GEFS_20261116_frio", "M4_SS__GEFS_20261116_central", "M1_SS__GEFS_20261116_central"):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True)
    for T in Ts:
        s = json.loads((SAL / "cache" / f"{CASOS[T][1]}__{e}.json").read_text(encoding="utf-8"))["serie"]; t = s["tiempo_h"]
        ax[0].plot(t, s["temperatura_maxima_c"], color=cols[T], label=f"{T} °C"); ax[1].plot(t, s["temperatura_minima_S25_c"], color=cols[T]); ax[2].plot(t, s["diferencia_maxima_S25_c"], color=cols[T])
    ax[0].axhline(70, color="k"); ax[2].axhline(20, color="k"); ax[0].legend()
    ax[0].set_title("Tmax [°C] (límite 70)"); ax[1].set_title("Tmin S25 [°C]"); ax[2].set_title("ΔT [°C] (límite 20)")
    for a in ax: a.set_xlabel("Tiempo desde t0 [h]"); a.grid(alpha=.3); a.set_xlim(0, 336)
    ax[0].set_ylabel("Temperatura [°C]"); ax[1].set_ylabel("Temperatura [°C]"); ax[2].set_ylabel("Diferencia térmica ΔT [°C]")
    fig.suptitle(f"Plan 1 con T de colocación 10–13 °C — escenario {e}"); fig.tight_layout()
    fig.savefig(SAL / "figuras" / f"curvas_T10_13_{e}.png", dpi=130); plt.close(fig)
