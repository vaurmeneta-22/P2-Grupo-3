"""Grafica las ocho curvas de Qdot del archivo de términos fuente."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/hidratacion/source_terms.txt"
OUTPUT = ROOT / "results/dossier_final/curvas_hidratacion_8_series.png"


def main() -> None:
    datos = np.loadtxt(SOURCE, comments="#")
    t = datos[:, 0]
    colores = {"M1": "#1769aa", "M2": "#e07a1f", "M3": "#2a8f62", "M4": "#8c5fb3"}
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    for mix, color in colores.items():
        indice = int(mix[1:])
        ax.plot(t, datos[:, indice], color=color, lw=2.2, label=f"{mix} · SS")
        ax.plot(t, datos[:, indice + 4], color=color, lw=1.8, ls="--", label=f"{mix} · 22A")
    ax.set(title="Curvas de generación volumétrica de calor", xlabel="Edad desde colocación [h]",
           ylabel="Q̇ [W/m³]")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, frameon=False)
    ax.text(0.01, 0.98, "Línea continua: SS · línea discontinua: 22A",
            transform=ax.transAxes, va="top", ha="left", fontsize=9, color="#4b5563")
    fig.savefig(OUTPUT, dpi=180, facecolor="white")
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
