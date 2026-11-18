"""Matriz de ΔT máximo (336 h): 11 casos × 24 escenarios = 264 valores.

Lee solo tablas ya calculadas de reports/Defensa/bateria7/tablas y escribe
``matriz_delta_11_casos.csv`` (formato ancho), ``matriz_delta_11_casos_larga.csv`` y la figura
``figuras/matriz_delta_11_casos.png``. No ejecuta simulaciones.
"""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

RAIZ = Path(__file__).resolve().parents[1]
TAB = RAIZ / "reports/Defensa/bateria7/tablas"
FIG = RAIZ / "reports/Defensa/bateria7/figuras"
CURVAS = ["M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"]
CLIMAS = [("frio", "Clima frío"), ("central", "Clima central"), ("calido", "Clima cálido")]
PESO_CLIMA = {"frio": 0.354838709677, "central": 0.322580645161, "calido": 0.322580645161}
CASOS = [
    (1, "Base original IFC, 18 °C", "base"), (2, "Base simultáneo, 18 °C", "base"), (3, "Base largo 24 h, 18 °C", "base"),
    (11, "Base original IFC a 10 °C", "base"),
    (7, "Control: 10 h sin manta, 10 °C", "control"),
    (4, "Plan 1: 10 h, PUR 60 mm, 10 °C", "plan"), (5, "Plan simultáneo, PUR 60 mm, 10 °C", "plan"), (6, "Plan largo 24 h, PUR 60 mm, 10 °C", "plan"),
    (8, "Plan 1 a 11 °C", "temp"), (9, "Plan 1 a 12 °C", "temp"), (10, "Plan 1 a 13 °C", "temp"),
]


def esc(c, w): return f"{c}__GEFS_20261116_{w}"


val = {}
nombres7 = {"1 Base original": 1, "2 Base simultáneo": 2, "3 Base largo 24 h": 3, "4 Plan 1 original 10 h": 4, "5 Plan simultáneo": 5, "6 Plan largo 24 h": 6, "7 Control sin manta": 7}
for r in csv.DictReader((TAB / "resultados_por_escenario.csv").open(encoding="utf-8")):
    val[(nombres7[r["caso"]], r["escenario"])] = float(r["M_delta_c"])
for r in csv.DictReader((TAB / "plan1_T10_13_por_escenario.csv").open(encoding="utf-8")):
    t = int(r["T_colocacion_c"])
    if t in (11, 12, 13):
        val[({11: 8, 12: 9, 13: 10}[t], r["escenario"])] = float(r["M_delta_c"])
for r in csv.DictReader((TAB / "caso11_por_escenario.csv").open(encoding="utf-8")):
    val[(11, r["escenario"])] = float(r["MD_caso11"])
assert len(val) == 264, len(val)
# consistencia: el caso 4 en plan1_T10_13 (T=10) y el de resultados_por_escenario coinciden
for r in csv.DictReader((TAB / "plan1_T10_13_por_escenario.csv").open(encoding="utf-8")):
    if int(r["T_colocacion_c"]) == 10:
        assert abs(val[(4, r["escenario"])] - float(r["M_delta_c"])) < 1e-3

cols = [(c, w) for w, _ in CLIMAS for c in CURVAS]
with (TAB / "matriz_delta_11_casos.csv").open("w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["caso", "nombre"] + [esc(c, w) for c, w in cols] + ["delta_max_c", "n_excedencias_de_24", "Pf_conjunta"])
    for n, nom, _ in CASOS:
        fila = [val[(n, esc(c, w))] for c, w in cols]
        pf = sum(PESO_CLIMA[w] / 8 for (c, w), v in zip(cols, fila) if v > 20)
        wr.writerow([n, nom] + [f"{v:.3f}" for v in fila] + [f"{max(fila):.3f}", sum(v > 20 for v in fila), f"{pf:.6f}"])
with (TAB / "matriz_delta_11_casos_larga.csv").open("w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(["caso", "nombre", "curva", "clima", "delta_max_336h_c", "excede_20"])
    for n, nom, _ in CASOS:
        for c, w in cols:
            v = val[(n, esc(c, w))]; wr.writerow([n, nom, c, w, f"{v:.3f}", int(v > 20)])
json.dump({"casos": [(n, nom, g) for n, nom, g in CASOS], "cols": cols,
           "val": {f"{n}|{c}|{w}": val[(n, esc(c, w))] for n, _, _ in CASOS for c, w in cols}},
          (TAB / "matriz_delta_11_casos.json").open("w", encoding="utf-8"))

M = np.array([[val[(n, esc(c, w))] for c, w in cols] for n, _, _ in CASOS])
from matplotlib.colors import TwoSlopeNorm
vmin, vmax = float(np.floor(M.min())), float(np.ceil(M.max()))
norm = TwoSlopeNorm(vmin=vmin, vcenter=20.0, vmax=vmax)
cmap = plt.get_cmap("RdYlBu_r")
fig, ax = plt.subplots(figsize=(21, 8.2))
im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]; r, g, b, _ = cmap(norm(v)); lum = 0.299 * r + 0.587 * g + 0.114 * b
        ax.text(j, i, f"{v:.1f}" + (" ×" if v > 20 else ""), ha="center", va="center", fontsize=9, color="black" if lum > 0.5 else "white", fontweight="bold" if v > 20 else None)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c for c, _ in cols], rotation=60, ha="right", fontsize=9)
ax.set_yticks(range(len(CASOS))); ax.set_yticklabels([f"{n} · {nom}" for n, nom, _ in CASOS], fontsize=10)
for k in (8, 16): ax.axvline(k - 0.5, color="white", lw=4)
for k in (4, 5, 8): ax.axhline(k - 0.5, color="white", lw=3)
for k, (w, lab) in enumerate(CLIMAS): ax.text(k * 8 + 3.5, -0.75, lab, ha="center", fontsize=12, fontweight="bold")
cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.022); cb.set_label("ΔT máximo en 336 h [°C] (azul < 20 °C · rojo > 20 °C)")
fig.suptitle("ΔT máximo (Tmax − Tmin S25) en 336 h por caso y escenario · 11 casos × 24 escenarios = 264 valores · × = supera 20 °C", fontsize=12, y=0.985)
fig.subplots_adjust(left=0.2, right=0.93, top=0.88, bottom=0.14); fig.savefig(FIG / "matriz_delta_11_casos.png", dpi=110); plt.close(fig)
print(vmin, vmax, "excedencias:", int((M > 20).sum()), "de", M.size)
