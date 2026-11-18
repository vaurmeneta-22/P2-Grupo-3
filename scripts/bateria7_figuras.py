"""Figuras de la batería de 7 casos: curvas comparables y vistas 3D con escala común."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
CASOS = [("bat1_base_original", "1 Base original", "#7f7f7f", "--"), ("bat2_base_simultaneo", "2 Base simultáneo", "#bcbd22", "--"),
         ("bat3_base_largo_24h", "3 Base largo 24 h", "#8c564b", "--"), ("bat4_plan1_original_10h", "4 Plan 1 (10 h)", "#1f77b4", "-"),
         ("bat5_plan_simultaneo", "5 Plan simultáneo", "#2ca02c", "-"), ("bat6_plan_largo_24h", "6 Plan largo 24 h", "#d62728", "-"),
         ("bat7_control_sin_manta", "7 Control sin manta", "#9467bd", ":")]
ESCS = ["M1_SS__GEFS_20261116_central", "M4_SS__GEFS_20261116_frio", "M4_SS__GEFS_20261116_central"]


def serie(p, e):
    return json.loads((SAL / "cache" / f"{p}__{e}.json").read_text(encoding="utf-8"))["serie"]


for e in ESCS:
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), sharex=True)
    for p, n, c, ls in CASOS:
        s = serie(p, e); t = np.asarray(s["tiempo_h"])
        ax[0].plot(t, s["temperatura_maxima_c"], color=c, ls=ls, label=n)
        ax[1].plot(t, s["temperatura_minima_S25_c"], color=c, ls=ls)
        ax[2].plot(t, s["diferencia_maxima_S25_c"], color=c, ls=ls)
    ax[0].axhline(70, color="k", lw=1); ax[2].axhline(20, color="k", lw=1)
    ax[0].set_title("Tmax del hormigón activo [°C] (límite 70)"); ax[1].set_title("Tmin S25 [°C]"); ax[2].set_title("ΔT = Tmax − Tmin S25 [°C] (límite 20)")
    for a in ax: a.set_xlabel("horas desde t0"); a.grid(alpha=.3); a.set_xlim(0, 336)
    ax[0].legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Escenario {e} (ΔT instantáneo; el máximo del criterio es el máximo de la curva)")
    fig.tight_layout(); fig.savefig(SAL / "figuras" / f"curvas_{e}.png", dpi=130); plt.close(fig)


def vista(e, t_obj, nombre):
    fig = plt.figure(figsize=(18, 9)); vmin, vmax = 10, 60
    for k, (p, n, *_ ) in enumerate(CASOS):
        d = np.load(SAL / "campos" / f"{p}__{e}.npz")
        i = int(np.argmin(np.abs(d["t_h"] - t_obj)))
        m = d["mascara"]; idx = np.argwhere(m); T = d["T"][i]; paso = d["paso"]; org = d["origen"]
        xyz = org + (idx + 0.5) * paso
        activo = np.isfinite(T)
        # corte: se retira la mitad y>ymed para ver el interior
        ym = 0.5 * (xyz[:, 1].min() + xyz[:, 1].max())
        sel = xyz[:, 1] <= ym
        ax = fig.add_subplot(2, 4, k + 1, projection="3d")
        sc = ax.scatter(xyz[sel, 0], xyz[sel, 1], xyz[sel, 2], c=T[sel], cmap="inferno", vmin=vmin, vmax=vmax, s=5, marker="s", linewidths=0)
        ax.set_title(f"{n}\nt={d['t_h'][i]:.0f} h  máx {np.nanmax(T):.1f} °C", fontsize=9)
        ax.view_init(24, 125); ax.set_box_aspect((np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]) / 2 + 1e-3, np.ptp(xyz[:, 2])))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    cax = fig.add_axes([0.62, 0.15, 0.3, 0.025]); fig.colorbar(sc, cax=cax, orientation="horizontal", label="T [°C] (escala común 10–60)")
    fig.suptitle(f"{e}: corte longitudinal (se omite la mitad y>mitad)"); fig.savefig(SAL / "figuras" / nombre, dpi=120); plt.close(fig)


vista("M4_SS__GEFS_20261116_frio", 164, "3d_M4_SS_frio_164h.png")
vista("M4_SS__GEFS_20261116_frio", 336, "3d_M4_SS_frio_336h.png")
vista("M4_SS__GEFS_20261116_central", 210, "3d_M4_SS_central_210h.png")
print("ok")
