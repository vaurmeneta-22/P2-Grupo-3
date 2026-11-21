"""Grafica el diferencial de los protocolos de retiro para M4_SS central."""
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt


RAIZ = Path(__file__).resolve().parents[1]
ENTRADA = RAIZ / "results/cierre_plan1"
SALIDA = RAIZ / "reports/Defensa/cierre_plan1/retiro_pur_100_120.png"
ESCENARIO = "M4_SS__GEFS_20261116_central.json.gz"


def main():
    figura, ejes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, sharey=True,
                               layout="constrained")
    for espesor, eje in zip((100, 120), ejes):
        for caso, titulo, color in (
            ("mantenido_1008", "Manta mantenida", "#246a73"),
            ("total_672", "Retiro total a 672 h", "#a53e3e"),
            ("capas_672", "Retiro en cuatro capas", "#a26910"),
        ):
            archivo = ENTRADA / f"pur_{espesor}mm_{caso}" / ESCENARIO
            with gzip.open(archivo, "rt", encoding="utf-8") as datos:
                serie = json.load(datos)["serie"]
            eje.plot([r["tiempo_h"] for r in serie],
                     [r["diferencia_maxima_S25_c"] for r in serie],
                     color=color, label=titulo)
        eje.axhline(20, color="#46505d", lw=1, ls="--", label="Umbral de 20 °C")
        eje.axvline(733, color="#46505d", lw=1, ls=":", label="Fin común de GEFS")
        eje.set_xlim(600, 1008)
        eje.set_ylim(17.5, 35.5)
        eje.set_title(f"PUR de {espesor} mm", loc="left")
        eje.set_ylabel("Diferencial S25 [°C]")
        eje.grid(axis="y", alpha=.2)
    ejes[0].legend(ncol=2, loc="upper right", fontsize=8)
    ejes[-1].set_xlabel("Horas desde el inicio del Plan 1")
    figura.suptitle("Retiro de PUR · escenario M4_SS con clima central")
    figura.savefig(SALIDA, dpi=160)
    plt.close(figura)


if __name__ == "__main__":
    main()
