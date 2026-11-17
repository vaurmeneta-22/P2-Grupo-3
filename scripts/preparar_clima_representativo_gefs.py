"""Prepara una única trayectoria representativa del pronóstico GEFS guardado.

Genera un snapshot climático completo, una configuración para el evaluador y
reemplaza la lista de climas del dossier por esa trayectoria con peso 1.
Desde la raíz del repo: python .\\scripts\\preparar_clima_representativo_gefs.py
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import demo_pipeline_ifc_3d as pipe
from mcplan.weather.serie import a_utc, iso_utc

DOSSIER = ROOT / "configs/evaluacion_final.jsonc"


def leer_jsonc(path: Path) -> dict:
    return json.loads(pipe.quitar_comentarios_jsonc(path.read_text(encoding="utf-8")))


def main(ensemble_path: Path | None = None) -> None:
    ensemble = Path(ensemble_path) if ensemble_path else (
        ROOT / "data/weather/2026-11-16T213220Z_GEFS_ensemble")
    if not ensemble.is_absolute():
        ensemble = ROOT / ensemble
    datos = json.loads((ensemble / "original.json").read_text(encoding="utf-8"))
    hourly = datos["hourly"]
    times = [a_utc(value) for value in hourly["time"]]
    keys = ["temperature_2m"] + [f"temperature_2m_member{i:02d}" for i in range(1, 31)]
    members = np.asarray([
        [np.nan if value is None else float(value) for value in hourly[key]]
        for key in keys
    ], dtype=float)
    counts = np.isfinite(members).sum(axis=0)

    dossier = leer_jsonc(DOSSIER)
    start = a_utc(dossier["inicio_obra_utc"])
    end = start + timedelta(hours=float(dossier["tiempo_final_h"]))
    target = [i for i, t in enumerate(times) if start <= t <= end]
    expected = int(round((end - start).total_seconds() / 3600)) + 1
    if len(target) != expected or times[target[0]] != start or times[target[-1]] != end:
        raise ValueError("El pronóstico no contiene todas las horas del período de simulación")
    if np.any(counts[target] != len(keys)):
        raise ValueError("Hay miembros sin datos dentro del período requerido para la simulación")

    # Para el gráfico se conserva todo el intervalo continuo común; el solver
    # utiliza únicamente la ventana t0–t_final configurada arriba.
    first = next(i for i, count in enumerate(counts) if count == len(keys))
    last = first
    while last + 1 < len(counts) and counts[last + 1] == len(keys):
        last += 1
    common_slice = slice(first, last + 1)
    curve = np.median(members[:, common_slice], axis=0)
    snapshot_times = times[common_slice]

    raw = (ensemble / "original.json").read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    ensemble_info = leer_jsonc(ensemble / "ensemble_metadata.json")
    issued = datetime.fromisoformat(ensemble_info["retrieved_at"].replace("Z", "+00:00"))
    date_id = issued.strftime("%Y%m%d")
    climate_id = f"GEFS_{date_id}_representativo"
    config_rel = f"configs/climas/gefs_{date_id}/representativo.jsonc"
    destination_name = ensemble.name.replace("_GEFS_ensemble", "_GEFS_representativo")
    destination = ROOT / "data/weather" / destination_name
    destination.mkdir(exist_ok=True)
    forecast = {
        "time_utc": [iso_utc(t) for t in snapshot_times],
        "temperature_2m": [round(float(v), 4) for v in curve],
        "units": {"temperature_2m": "degC"},
    }
    (destination / "forecast.json").write_text(
        json.dumps(forecast, indent=2), encoding="utf-8")
    metadata = {
        "source": "Open-Meteo Ensemble API; NOAA Global Ensemble Forecast System",
        "model": "NOAA GEFS 0.5° (ncep_gefs05)",
        "member": "representative_pointwise_median",
        "derivation": "Median temperature at each timestamp across the complete ensemble.",
        "ensemble_member_count": len(keys),
        "site": "Grupo 3: Mantos de Hornillos / Ovalle, Camarico – Ecos del Mar, Coquimbo",
        "latitude": -31.1115,
        "longitude": -71.6147,
        "timezone": "UTC",
        "coordinates_source": "Tabla docente de sitios del curso; mismos datos usados en E11",
        "issued_at": None,
        "issued_at_note": "la API no expone la hora de emisión del ciclo GEFS",
        "retrieved_at": ensemble_info["retrieved_at"],
        "valid_from": iso_utc(snapshot_times[0]),
        "valid_to": iso_utc(snapshot_times[-1]),
        "variables": ["temperature_2m"],
        "units": {"temperature_2m": "degC"},
        "original_file": f"../{ensemble.name}/original.json",
        "original_sha256": digest,
        "n_points": len(snapshot_times),
        "series_kind": "pronostico",
        "synthetic": False,
        "native_temporal_resolution": "3-hourly GEFS output; hourly time series interpolated by Open-Meteo",
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    climate_config = {
        "clima": {
            "fuente": "snapshot",
            "inicio_obra_utc": iso_utc(start),
            "snapshot": {
                "carpeta": destination.relative_to(ROOT).as_posix(),
                "raiz": "data/weather",
                "extension": None,
            },
        }
    }
    config_path = ROOT / config_rel
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(climate_config, indent=2) + "\n", encoding="utf-8")

    block = (
        '"climas": [\n'
        '    // Una trayectoria climática de referencia; peso total 1.\n'
        f'    {{"id": "{climate_id}", "config_jsonc": "{config_rel}", "peso": 1.0}}\n'
        '  ]'
    )
    text = DOSSIER.read_text(encoding="utf-8")
    updated, count = re.subn(
        r'"climas": \[[\s\S]*?\n  \],(?=\n  "pesos_curva_clima_independientes")',
        block + ",", text, count=1)
    if count != 1:
        raise ValueError("No se localizó la lista climas del dossier")
    updated = re.sub(
        r'"escenario_id":\s*"[^"]+"',
        f'"escenario_id": "M1_SS__{climate_id}"', updated, count=1)
    DOSSIER.write_text(updated, encoding="utf-8")

    print(f"Snapshot climático preparado: {destination.relative_to(ROOT)}")
    print(f"Cobertura del snapshot: {iso_utc(snapshot_times[0])} a {iso_utc(snapshot_times[-1])}")
    print(f"Ventana de simulación: {iso_utc(start)} a {iso_utc(end)} ({expected} datos)")
    print(f"Configuración del dossier: {climate_id}, peso 1.0")


if __name__ == "__main__":
    try:
        main(Path(sys.argv[1]) if len(sys.argv) > 1 else None)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
