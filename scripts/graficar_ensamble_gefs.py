"""Grafica la trayectoria de temperatura del snapshot de clima usado por el dossier.

Desde la raíz del repositorio:
    python .\\scripts\\graficar_ensamble_gefs.py
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "data/weather/2026-11-16T213220Z_GEFS_representativo"
DOSSIER = ROOT / "configs/evaluacion_final.jsonc"
OUTPUT = ROOT / "results/dossier_final/pronostico_gefs_referencia_14_dias.png"
UTC = timezone.utc


def main() -> None:
    forecast = json.loads((SNAPSHOT / "forecast.json").read_text(encoding="utf-8"))
    metadata = json.loads((SNAPSHOT / "metadata.json").read_text(encoding="utf-8"))
    times = [datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(UTC)
             for value in forecast["time_utc"]]
    temperature = forecast["temperature_2m"]
    if len(times) != len(temperature) or len(times) < 2:
        raise ValueError("El snapshot no contiene una serie horaria válida")

    import demo_pipeline_ifc_3d as pipe
    sys.path.insert(0, str(ROOT / "scripts"))
    dossier = json.loads(pipe.quitar_comentarios_jsonc(
        DOSSIER.read_text(encoding="utf-8")))
    start = datetime.fromisoformat(dossier["inicio_obra_utc"].replace("Z", "+00:00"))
    end = start + timedelta(hours=float(dossier["tiempo_final_h"]))
    if times[0] > start or times[-1] < end:
        raise ValueError("El snapshot no cubre toda la ventana simulada del dossier")

    fig, ax = plt.subplots(figsize=(12, 5.8), layout="constrained")
    ax.plot(times, temperature, color="#1976c9", linewidth=1.8,
            label="Temperatura del aire del snapshot", zorder=3)
    ax.axvspan(start, end, color="#ffc857", alpha=0.28,
               label=f"Ventana simulada ({dossier['tiempo_final_h']:g} h desde el inicio de obra)",
               zorder=1)
    ax.axvline(times[0], color="#555555", linestyle=":", linewidth=1.5, zorder=2)
    ax.axvline(times[-1], color="#555555", linestyle=":", linewidth=1.5, zorder=2)

    coverage_handle = Line2D([0], [0], color="#555555", linestyle=":", linewidth=1.5,
                             label="Intervalo cubierto por el snapshot")
    handles, labels = ax.get_legend_handles_labels()
    handles.append(coverage_handle)
    labels.append("Intervalo cubierto por el snapshot")
    ax.legend(handles, labels, loc="lower right", framealpha=0.94, fontsize=9)

    issued = datetime.fromisoformat(metadata["retrieved_at"].replace("Z", "+00:00"))
    snapshot_label = issued.strftime("%Y-%m-%dT%H%M%SZ")
    ax.set_title(f"Aire (snapshot {snapshot_label}, tipo: pronóstico)", fontsize=14, pad=12)
    ax.set_xlabel("Fecha y hora (UTC)")
    ax.set_ylabel("T aire [°C]")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2, tz=UTC))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m %H:%M", tz=UTC))
    ax.grid(True, color="#d8d8d8", alpha=0.6)
    ax.tick_params(axis="x", labelrotation=32)
    ax.margins(x=0.015)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180, facecolor="white")
    plt.close(fig)
    print(f"Gráfico guardado: {OUTPUT}")
    print(f"Ventana simulada: {start.isoformat()} a {end.isoformat()}")
    print(f"Temperatura en la trayectoria: {min(temperature):.1f} a {max(temperature):.1f} °C")


if __name__ == "__main__":
    main()
