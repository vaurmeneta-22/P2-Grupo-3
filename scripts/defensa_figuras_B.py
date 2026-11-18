"""Figuras de la propuesta B y del retiro cada 14 días (falla) a partir de results/defensa/cache2.

    python scripts/defensa_figuras_B.py          # escribe reports/Defensa/figuras/propuestaB_*.png
"""
import glob, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
FIG = RAIZ / "reports/Defensa/figuras"
CACHE2 = RAIZ / "results/defensa/cache2"
VOL = {1: 147.626, 2: 147.626, 3: 214.224, 4: 88.925, 5: 38.375}
DIAS_B, DIAS_A = (56, 77, 98, 119), (56, 70, 84, 98)


def cargar(patron):
    f = glob.glob(str(CACHE2 / patron))
    if not f:
        raise FileNotFoundError(patron)
    return json.load(open(f[0], encoding="utf-8"))


def figura_series():
    fig, ax = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
    for col, (suelo, suf) in enumerate((("suelo efectivo 2 m", "L2m"), ("suelo efectivo 0,8 m", "L0_8m"))):
        B = cargar(f"retB_{suf}__M4_SS__calida__*.json")
        A = cargar(f"ret2_A_14d_desde56_(falla)_{suf}__M4_SS__calida__*.json")
        for d, etq, c, dias in ((A, "retiro cada 14 días (falla)", "#c44e52", DIAS_A), (B, "propuesta B: retiro cada 21 días", "#4c72b0", DIAS_B)):
            s = d["serie"]; t = np.asarray(s["h"]) / 24
            ax[0][col].plot(t, s["delta_c"], color=c, label=etq)
            ax[1][col].plot(t, s["Tmax_c"], color=c, label=etq)
            for x in dias:
                ax[0][col].axvline(x, color=c, ls=":", lw=0.8)
        ax[0][col].axhline(20, color="k", ls="--"); ax[0][col].text(1, 20.3, "20 °C", fontsize=8)
        ax[1][col].axhline(70, color="k", ls="--"); ax[1][col].text(1, 71, "70 °C", fontsize=8)
        ax[0][col].set(title=f"M4_SS, clima cálido — {suelo}", ylabel="ΔT = Tmax − T_S25,min [°C]", xlim=(0, 141))
        ax[1][col].set(xlabel="días desde t0 (19-11-2026 10:00, Ovalle)", ylabel="Tmax del hormigón [°C]", ylim=(0, 75))
        ax[0][col].legend(fontsize=8); ax[0][col].grid(alpha=0.3); ax[1][col].grid(alpha=0.3)
    fig.suptitle("Candidato 36 m³/h progresivo, 10 °C, PUR 80 mm, aislación inferior 100 mm — extensión climática supuesta desde el día 28\n"
                 "(líneas punteadas: días de retiro de cada capa de 20 mm)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "propuestaB_series_B_vs_retiro14d.png", dpi=140); plt.close(fig)


def figura_cronograma():
    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5))
    t = 0.0; col = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974"]
    for i in range(1, 6):
        dur = VOL[i] / 36.0
        ax[0].barh(0, dur, left=t, color=col[i - 1], edgecolor="k")
        ax[0].text(t + dur / 2, 0, f"V{i}\n{dur:.1f} h", ha="center", va="center", fontsize=7, color="w")
        t += dur
    ax[0].axvline(15, color="darkorange"); ax[0].text(15.2, 0.35, "manta 80 mm desde las 15 h", fontsize=8, color="darkorange")
    ax[0].set(xlim=(0, 24), yticks=[], xlabel="horas desde t0", title=f"Colocación continua progresiva a 36 m³/h (una cuadrilla, sin esperas): D = {t:.2f} h")
    esp = [(0, 56, 80), (56, 77, 60), (77, 98, 40), (98, 119, 20), (119, 141, 0)]
    for a, b, e in esp:
        ax[1].barh(0, b - a, left=a, color=plt.cm.Blues(0.15 + e / 100), edgecolor="k")
        ax[1].text((a + b) / 2, 0, f"{e} mm" if e else "sin manta", ha="center", va="center", fontsize=8)
    for d in DIAS_B:
        ax[1].axvline(d, color="crimson", ls="--"); ax[1].text(d + 0.5, 0.42, f"−20 mm\nd{d}", fontsize=7, color="crimson")
    ax[1].set(xlim=(0, 141), yticks=[], xlabel="días desde t0 (19-11-2026 10:00)",
              title="Espesor de manta PUR sobre caras laterales y superior (aislación inferior de 100 mm permanente; moldaje 18 mm en V1, V2, V5 hasta 72 h)")
    fig.tight_layout(); fig.savefig(FIG / "propuestaB_cronograma.png", dpi=140); plt.close(fig)


if __name__ == "__main__":
    figura_series(); figura_cronograma(); print("ok")
