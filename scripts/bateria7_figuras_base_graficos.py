"""Gráficos grandes de «Caso base»: por variante, (a) Tmax y Tmin S25 superpuestas y (b) ΔT.

Mismo escenario (M4_SS / clima frío) y mismos tiempos (0–336 h) en ambos. Solo grafica series ya
calculadas (results/bateria7/cache); no ejecuta simulaciones.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
OUT = RAIZ / "reports/Defensa/bateria7/figuras"
ESC = "M4_SS__GEFS_20261116_frio"
VAR = [("v1", "bat1_base_original", "Original IFC · 18 °C"), ("v2", "bat2_base_simultaneo", "Simultáneo (t = 0) · 18 °C"),
       ("v3", "bat3_base_largo_24h", "Largo 24 h · 18 °C"), ("v4", "bat11_base_original_T10", "Original IFC a 10 °C")]
C_TMAX, C_TMIN, C_DELTA = "#d7191c", "#2b83ba", "#6a3d9a"
plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 12})
ESCENARIO_TXT = "escenario M4_SS / clima frío"

for k, plan, nombre in VAR:
    s = json.loads((SAL / "cache" / f"{plan}__{ESC}.json").read_text(encoding="utf-8"))["serie"]
    t = np.asarray(s["tiempo_h"]); tmax = np.asarray(s["temperatura_maxima_c"]); tmin = np.asarray(s["temperatura_minima_S25_c"]); d = np.asarray(s["diferencia_maxima_S25_c"])
    # (a) temperaturas
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.plot(t, tmax, color=C_TMAX, lw=2.6, label="Tmax del hormigón")
    ax.plot(t, tmin, color=C_TMIN, lw=2.6, label="Tmin S25 (25 mm bajo la superficie)")
    ax.axhline(70, color="k", ls="--", lw=1.8, label="Referencia 70 °C")
    ax.set_xlim(0, 336); ax.set_ylim(0, 76); ax.set_xticks(range(0, 337, 48)); ax.grid(alpha=.3)
    ax.set_xlabel("Tiempo desde t0 [h]"); ax.set_ylabel("Temperatura [°C]")
    ax.set_title(f"{nombre}\nTmax y Tmin S25 · {ESCENARIO_TXT}", fontsize=13)
    ax.legend(loc="center right", bbox_to_anchor=(0.995, 0.47), framealpha=.92)
    fig.tight_layout(pad=0.6); fig.savefig(OUT / f"base_{k}_temperaturas.png", dpi=130); plt.close(fig)
    # (b) delta
    im = int(np.argmax(d))
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.plot(t, d, color=C_DELTA, lw=2.6, label="ΔT = Tmax − Tmin S25")
    ax.axhline(20, color="#e6550d", lw=3.2, label="Límite 20 °C")
    ax.plot(t[im], d[im], "o", color="k", ms=10, zorder=5, label=f"Máximo: {d[im]:.2f} °C a las {t[im]:.0f} h".replace(".", ","))
    ax.annotate(f"{d[im]:.2f} °C\na las {t[im]:.0f} h".replace(".", ","), (t[im], d[im]), xytext=(18, 14) if t[im] < 250 else (-110, 14), textcoords="offset points", fontsize=12, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color="k", lw=1))
    ax.set_xlim(0, 336); ax.set_ylim(0, 52); ax.set_xticks(range(0, 337, 48)); ax.grid(alpha=.3)
    ax.set_xlabel("Tiempo desde t0 [h]"); ax.set_ylabel("Diferencia térmica ΔT [°C]")
    ax.set_title(f"{nombre}\nΔT = Tmax − Tmin S25 · {ESCENARIO_TXT}", fontsize=13)
    ax.legend(loc="lower right", framealpha=.92)
    fig.tight_layout(pad=0.6); fig.savefig(OUT / f"base_{k}_delta.png", dpi=130); plt.close(fig)
    print(k, f"tmax_max_escenario={tmax.max():.2f}", f"delta_max={d[im]:.3f} a {t[im]:.1f} h")
