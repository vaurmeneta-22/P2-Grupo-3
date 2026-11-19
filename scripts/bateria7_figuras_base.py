"""Figuras de la pestaña «Caso base» de la presentación: 4 variantes sin manta (casos 1, 2, 3 y 11).

Solo grafica resultados ya calculados (results/bateria7/cache y campos); no ejecuta simulaciones.
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
CFG = json.loads((RAIZ / "configs/defensa/bateria7_casos.json").read_text(encoding="utf-8"))
PLANES = {p["id"]: p for p in CFG["planes"]}
VAR = [("v1", "bat1_base_original", "1 · Original IFC, 18 °C", "#178A57"), ("v2", "bat2_base_simultaneo", "2 · Simultáneo (t=0), 18 °C", "#D94F4F"),
       ("v3", "bat3_base_largo_24h", "3 · Largo 24 h, 18 °C", "#2678C8"), ("v4", "bat11_base_original_T10", "4 · Original IFC a 10 °C", "#B69500")]
ESC = {"frio": "M4_SS__GEFS_20261116_frio", "central": "M1_SS__GEFS_20261116_central"}
PCOL = ["#7be0bd", "#79c8bb", "#ffc17f", "#a9a0ff", "#f38ba8"]
plt.rcParams.update({"font.size": 10})


def serie(plan, e):
    return json.loads((SAL / "cache" / f"{plan}__{e}.json").read_text(encoding="utf-8"))["serie"]


# 1) curvas por variante
for k, plan, nombre, col in VAR:
    s = serie(plan, ESC["frio"]); t = s["tiempo_h"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.1), sharex=True)
    ax[0].plot(t, s["temperatura_maxima_c"], color=col, lw=2); ax[0].axhline(70, color="k"); ax[0].set_title("Tmax [°C] (límite 70)"); ax[0].set_ylim(0, 75)
    ax[1].plot(t, s["temperatura_minima_S25_c"], color=col, lw=2); ax[1].set_title("Tmin S25 [°C]")
    ax[2].plot(t, s["diferencia_maxima_S25_c"], color=col, lw=2); ax[2].axhline(20, color="k"); ax[2].set_title("ΔT = Tmax − Tmin S25 [°C] (límite 20)")
    imax = int(np.argmax(s["diferencia_maxima_S25_c"])); ax[2].plot(t[imax], s["diferencia_maxima_S25_c"][imax], "o", color="k")
    ax[2].annotate(f"máx {s['diferencia_maxima_S25_c'][imax]:.2f} °C a {t[imax]:.0f} h", (t[imax], s["diferencia_maxima_S25_c"][imax]), xytext=(-150, -28), textcoords="offset points", fontsize=9)
    for a in ax: a.set_xlabel("horas desde t0"); a.grid(alpha=.3); a.set_xlim(0, 336)
    fig.suptitle(f"Variante {nombre} — escenario M4_SS / clima frío (el de mayor ΔT)"); fig.tight_layout()
    fig.savefig(OUT / f"base_{k}_curvas.png", dpi=130); plt.close(fig)

# 2) superposición ΔT, un escenario por figura
for etiq, e in ESC.items():
    fig, ax = plt.subplots(figsize=(12, 5.5), facecolor="#fbfcfe")
    ax.set_facecolor("#fbfcfe")
    ax.axhspan(20, 50, color="#d9534f", alpha=.055, zorder=0)
    for k, plan, nombre, col in VAR:
        s = serie(plan, e)
        label = nombre.replace("Original IFC,", "Original IFC ·").replace("Simultáneo (t=0),", "Simultáneo ·").replace("Largo 24 h,", "Largo 24 h ·").replace("Original IFC a", "Original IFC ·")
        ax.plot(s["tiempo_h"], s["diferencia_maxima_S25_c"], color=col, lw=2.8, ls="--" if k == "v3" else "-", label=label, zorder=3)
    ax.axhline(20, color="#c44848", lw=1.8, ls=(0, (5, 4)), zorder=2)
    ax.text(4, 21.1, "Límite ΔT · 20 °C", color="#a83232", fontsize=9, weight="bold",
            bbox={"boxstyle": "round,pad=.35", "facecolor": "white", "edgecolor": "#f1caca", "alpha": .96})
    ax.set_xlabel("Tiempo desde el inicio de obra [h]", labelpad=9)
    ax.set_ylabel("Diferencia térmica · ΔT [°C]", labelpad=9)
    ax.set_xlim(0, 336); ax.set_ylim(0, 50); ax.set_xticks(np.arange(0, 337, 48))
    ax.grid(axis="both", color="#cbd5e1", alpha=.55, linewidth=.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#64748b"); ax.spines["bottom"].set_color("#64748b")
    ax.tick_params(colors="#334155", labelsize=9)
    escenario = e.replace("M4_SS__GEFS_20261116_frio", "M4_SS · clima frío").replace("M1_SS__GEFS_20261116_central", "M1_SS · clima central")
    ax.set_title(f"Evolución de ΔT · cuatro variantes sin manta\nEscenario {escenario}", color="#172033", fontsize=14, weight="bold", pad=14)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2, frameon=False,
              fontsize=9, labelcolor="#263445", columnspacing=2.0, handlelength=3.0)
    fig.subplots_adjust(left=.09, right=.985, top=.84, bottom=.25)
    fig.savefig(OUT / f"base_superposicion_delta_{etiq}.png", dpi=180, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)

# 3) cronogramas
for k, plan, nombre, col in VAR:
    p = PLANES[plan]; act = {int(i): v for i, v in p["activaciones_h"].items()}; dur = {int(i): v for i, v in p["duraciones_h"].items()}
    fig, ax = plt.subplots(figsize=(10, 3.4))
    for i in range(1, 6):
        simult = k == "v2"
        ax.barh(6 - i, dur[i], left=act[i], color=PCOL[i - 1], edgecolor="k", linewidth=.7, hatch="//" if simult else None, alpha=.55 if simult else 1)
        ax.plot(act[i], 6 - i, "v", color="k", ms=7)
        ax.text(act[i] + dur[i] + .3, 6 - i, f"V{i}: inicio {act[i]:.2f} h" + (f", dur. {dur[i]:.2f} h" if not simult else ""), va="center", fontsize=8)
        if i in (1, 2, 5):
            ax.plot(act[i] + 72, 6 - i, "s", color="#c0392b", ms=6)
    ax.set_yticks(range(1, 6)); ax.set_yticklabels([f"V{i}" for i in (5, 4, 3, 2, 1)]); ax.set_xlim(0, 120); ax.set_xlabel("horas desde t0 (19-11-2026 10:00)"); ax.grid(axis="x", alpha=.3)
    sub = ("▼ activación térmica instantánea al inicio de cada vertida · ■ rojo: retiro del moldaje (72 h tras colocar; solo V1, V2 y V5) · sin manta"
           + ("\nbarras rayadas: duración por capacidad (6,368 h), cota teórica; la activación simultánea es una idealización térmica, no una colocación ejecutable" if k == "v2" else ""))
    ax.set_title(f"Cronograma — {nombre}", fontsize=11)
    fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.text(0.01, 0.01, sub, fontsize=7.5, va="bottom"); fig.savefig(OUT / f"base_{k}_cronograma.png", dpi=130); plt.close(fig)

# 4) 3D por variante, mismo escenario, mismos tiempos y escala
for k, plan, nombre, col in VAR:
    d = np.load(SAL / "campos" / f"{plan}__{ESC['frio']}.npz"); idx = np.argwhere(d["mascara"]); xyz = d["origen"] + (idx + 0.5) * d["paso"]
    sel = xyz[:, 1] <= 0.5 * (xyz[:, 1].min() + xyz[:, 1].max())
    fig = plt.figure(figsize=(13, 4.6))
    for j, tobj in enumerate((164, 336)):
        i = int(np.argmin(np.abs(d["t_h"] - tobj))); T = d["T"][i]
        ax = fig.add_subplot(1, 2, j + 1, projection="3d")
        sc = ax.scatter(xyz[sel, 0], xyz[sel, 1], xyz[sel, 2], c=T[sel], cmap="inferno", vmin=10, vmax=60, s=6, marker="s", linewidths=0)
        ax.set_title(f"t = {d['t_h'][i]:.0f} h · Tmax {np.nanmax(T):.1f} °C", fontsize=10); ax.view_init(24, 125)
        ax.set_box_aspect((np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]) / 2 + 1e-3, np.ptp(xyz[:, 2]))); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    cax = fig.add_axes([0.3, 0.13, 0.4, 0.03]); fig.colorbar(sc, cax=cax, orientation="horizontal", label="T [°C] — escala común 10–60 °C")
    fig.suptitle(f"{nombre} — corte longitudinal, escenario M4_SS / clima frío", fontsize=11); fig.savefig(OUT / f"base_{k}_3d.png", dpi=120); plt.close(fig)
print("ok")
