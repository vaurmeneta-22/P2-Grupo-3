"""Grafica las tres trayectorias GEFS completas usadas en el dossier."""
from datetime import datetime
import json
from pathlib import Path
from datetime import timezone, timedelta

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias"
OUTPUT = ROOT / "results/dossier_final/clima_gefs_tres_trayectorias_14d.png"
SERIES = (("fria", "Fría", "#2878b5"),
          ("central", "Central", "#6b7280"),
          ("calida", "Cálida", "#d97706"))


def main() -> None:
    fig, ax = plt.subplots(figsize=(11, 5.4), constrained_layout=True)
    # El intervalo de noviembre de 2026 está íntegramente en UTC−3 en Chile.
    # Usar el desplazamiento explícito evita depender del paquete tzdata en Windows.
    zona_local = timezone(timedelta(hours=-3))
    for folder, label, color in SERIES:
        raw = json.loads((BASE / folder / "forecast.json").read_text(encoding="utf-8"))
        times_utc = [datetime.fromisoformat(t.replace("Z", "+00:00")) for t in raw["time_utc"]]
        times = [t.astimezone(zona_local) for t in times_utc]
        values = raw["temperature_2m"]
        if len(times) != 337 or len(values) != 337:
            raise ValueError(f"{folder}: se esperaban 337 pares horario/temperatura")
        ax.plot(times, values, color=color, lw=2.0, label=label)
    ax.set(title="Pronóstico de temperatura del aire · 14 días",
           xlabel="Fecha y hora local (Ovalle)", ylabel="T aire [°C]")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2, tz=zona_local))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m %H:%M", tz=zona_local))
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=3, loc="upper right")
    fig.autofmt_xdate(rotation=25, ha="right")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=180, facecolor="white")
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
