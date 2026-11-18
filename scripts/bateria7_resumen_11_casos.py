"""Resumen de los 11 casos (cuatro paneles): D, Pf conjunta, ΔT máximo y Tmax máxima.

Todo se lee de tablas ya calculadas en reports/Defensa/bateria7/tablas y de la configuración de
los planes; no se transcribe ningún valor a mano ni se ejecutan simulaciones. Escribe
``tablas/resumen_11_casos.csv`` y ``figuras/resumen_11_casos.png``.
"""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RAIZ = Path(__file__).resolve().parents[1]
TAB = RAIZ / "reports/Defensa/bateria7/tablas"
FIG = RAIZ / "reports/Defensa/bateria7/figuras"
CFG = json.loads((RAIZ / "configs/defensa/bateria7_casos.json").read_text(encoding="utf-8"))
PLAN_ID = {1: "bat1_base_original", 2: "bat2_base_simultaneo", 3: "bat3_base_largo_24h", 4: "bat4_plan1_original_10h", 5: "bat5_plan_simultaneo",
           6: "bat6_plan_largo_24h", 7: "bat7_control_sin_manta", 8: "bat4_T11", 9: "bat4_T12", 10: "bat4_T13", 11: "bat11_base_original_T10"}
SIMULT = {2, 5}
ORDEN = [1, 2, 3, 11, 7, 4, 5, 6, 8, 9, 10]          # mismo orden que la matriz de ΔT
COLOR = {1: "#8a8f94", 2: "#a9a28a", 3: "#7d6b5a", 11: "#d9a066",          # caso base (18 °C y 10 °C)
         7: "#9467bd",                                                      # control sin manta
         4: "#1f77b4", 5: "#2ca02c", 6: "#17becf",                          # planes protegidos
         8: "#f2a83b", 9: "#e8743b", 10: "#c0392b"}                         # Plan 1 a 11, 12 y 13 °C
GRUPO = {1: "base", 2: "base", 3: "base", 11: "base", 7: "control", 4: "plan", 5: "plan", 6: "plan", 8: "temp", 9: "temp", 10: "temp"}


def esc(c, w): return f"{c}__GEFS_20261116_{w}"


# --- nombres, ΔT máximo y Pf desde la matriz final
nombre, dmax, pf, nexc = {}, {}, {}, {}
for r in csv.DictReader((TAB / "matriz_delta_11_casos.csv").open(encoding="utf-8")):
    n = int(r["caso"]); nombre[n] = r["nombre"]; dmax[n] = float(r["delta_max_c"]); pf[n] = float(r["Pf_conjunta"]); nexc[n] = int(r["n_excedencias_de_24"])
# --- Tmax máxima por caso desde las tablas por escenario
tmax = {}
nom7 = {"1 Base original": 1, "2 Base simultáneo": 2, "3 Base largo 24 h": 3, "4 Plan 1 original 10 h": 4, "5 Plan simultáneo": 5, "6 Plan largo 24 h": 6, "7 Control sin manta": 7}
for r in csv.DictReader((TAB / "resultados_por_escenario.csv").open(encoding="utf-8")):
    n = nom7[r["caso"]]; tmax[n] = max(tmax.get(n, 0), float(r["M_T_c"]))
for r in csv.DictReader((TAB / "plan1_T10_13_por_escenario.csv").open(encoding="utf-8")):
    t = int(r["T_colocacion_c"])
    if t in (11, 12, 13):
        n = {11: 8, 12: 9, 13: 10}[t]; tmax[n] = max(tmax.get(n, 0), float(r["M_T_c"]))
for r in csv.DictReader((TAB / "caso11_por_escenario.csv").open(encoding="utf-8")):
    tmax[11] = max(tmax.get(11, 0), float(r["MT_caso11"]))
# --- D desde la configuración (fin de colocación programado)
D = {}
planes = {p["id"]: p for p in CFG["planes"]}
for n, pid in PLAN_ID.items():
    p = planes[pid]
    D[n] = max(float(p["activaciones_h"][k]) + float(p["duraciones_h"][k]) for k in p["activaciones_h"])
assert set(tmax) == set(PLAN_ID) and set(dmax) == set(PLAN_ID) and set(pf) == set(PLAN_ID)

with (TAB / "resumen_11_casos.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["caso", "nombre", "D_h", "D_es_cota_por_capacidad_idealizada", "Pf_conjunta", "Pf_pct", "delta_max_entre_escenarios_c", "Tmax_max_entre_escenarios_c", "escenarios_con_delta_mayor_20"])
    for n in ORDEN:
        w.writerow([n, nombre[n], f"{D[n]:.5f}", int(n in SIMULT), f"{pf[n]:.6f}", f"{pf[n]*100:.4f}", f"{dmax[n]:.3f}", f"{tmax[n]:.3f}", nexc[n]])

plt.rcParams.update({"font.size": 12})
fig, axs = plt.subplots(1, 4, figsize=(16, 8.4), sharey=True, gridspec_kw={"width_ratios": [1.0, 1.0, 1.0, 1.0], "wspace": 0.08})
y = np.arange(len(ORDEN))[::-1]
etiquetas = [f"{n} · {nombre[n]}" for n in ORDEN]
paneles = [
    ("Duración constructiva D [h]", [D[n] for n in ORDEN], None, "Duración constructiva D", (0, 30)),
    ("Probabilidad conjunta de fallo Pf [%]", [pf[n] * 100 for n in ORDEN], 5, "Probabilidad conjunta\nde fallo Pf", (0, 125)),
    ("Diferencial térmico máximo ΔT [°C]", [dmax[n] for n in ORDEN], 20, "ΔT máximo\n(entre escenarios)", (0, 60)),
    ("Temperatura máxima Tmax [°C]", [tmax[n] for n in ORDEN], 70, "Tmax máxima\n(entre escenarios)", (0, 80)),
]
for ax, (xl, vals, lim, titulo, xlim) in zip(axs, paneles):
    for yi, n, v in zip(y, ORDEN, vals):
        hatch = "///" if (n in SIMULT and ax is axs[0]) else None
        ax.barh(yi, v, color=COLOR[n], edgecolor="#333" if hatch else COLOR[n], hatch=hatch, height=0.72, alpha=0.55 if hatch else 1)
        if ax is axs[0]:
            txt = f"{v:.3f}".replace(".", ",") + " (cota)" if n in SIMULT else f"{v:.1f}".replace(".", ",")
        elif ax is axs[1]:
            txt = "0" if v == 0 else (f"{v:.2f}".replace(".", ",") if v < 100 else "100")
            txt = ("12,5" if abs(v - 12.5) < 0.05 else txt)
        else:
            txt = f"{v:.1f}".replace(".", ",")
        ax.text(v + (xlim[1] * 0.012), yi, txt, va="center", fontsize=10.5, fontweight="bold" if (lim and v > lim) else None, color="#b2182b" if (lim and v > lim) else "#222", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=6)
    if lim is not None:
        ax.axvline(lim, color="#d7191c", ls="--", lw=2.2)
        ax.text(lim, len(ORDEN) - 0.35, f" límite {lim:g} {'%' if ax is axs[1] else '°C'}", color="#d7191c", fontsize=10.5, va="bottom", ha="left" if lim < xlim[1] * 0.8 else "right")
    ax.set_xlim(*xlim); ax.set_xlabel(xl); ax.set_title(titulo, fontsize=13, fontweight="bold", pad=24)
    ax.grid(axis="x", alpha=.3); ax.set_ylim(-0.6, len(ORDEN) - 0.4)
    for lin in (3.5 - 0, ):
        pass
    # separadores de grupo: base | control | planes | temperatura
    for yy in (len(ORDEN) - 4 - 0.5, len(ORDEN) - 5 - 0.5, len(ORDEN) - 8 - 0.5):
        ax.axhline(yy, color="#666", lw=0.8)
axs[0].set_yticks(y); axs[0].set_yticklabels(etiquetas, fontsize=11); axs[0].set_ylabel("Caso (11 casos · 24 escenarios cada uno)")
fig.legend(handles=[Patch(facecolor="#ccc", edgecolor="#333", hatch="///", label="Barra rayada = cota por capacidad (636,776 m³ / 100 m³/h = 6,368 h), asociada a una activación térmica idealizada; no es una duración ejecutable")], loc="lower center", bbox_to_anchor=(0.5, 0.0), fontsize=10.5, frameon=False)
fig.suptitle("Comparación de los 11 casos · 264 combinaciones evaluadas (8 curvas × 3 climas por caso, 336 h)", fontsize=14, fontweight="bold", y=0.985)
fig.subplots_adjust(left=0.205, right=0.995, top=0.835, bottom=0.115)
fig.savefig(FIG / "resumen_11_casos.png", dpi=125); plt.close(fig)
json.dump({"orden": ORDEN, "color": COLOR, "nombre": nombre, "grupo": GRUPO}, (TAB / "resumen_11_casos_colores.json").open("w", encoding="utf-8"), ensure_ascii=False)
for n in ORDEN:
    print(n, nombre[n], f"D={D[n]:.3f}", f"Pf={pf[n]*100:.2f}%", f"dT={dmax[n]:.3f}", f"Tmax={tmax[n]:.3f}")
