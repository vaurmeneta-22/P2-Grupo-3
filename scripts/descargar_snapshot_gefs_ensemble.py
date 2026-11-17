"""Descarga el ensamble NOAA GEFS y prepara su trayectoria climática representativa.

La respuesta cruda y los miembros individuales se archivan para trazabilidad. El dossier
usa una sola trayectoria obtenida como mediana horaria del ensamble.

Ejecutar desde la raíz:
    python .\\scripts\\descargar_snapshot_gefs_ensemble.py
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts"))
sys.path.insert(0, str(RAIZ / "src"))

import demo_pipeline_ifc_3d as pipe
from mcplan.weather.serie import a_utc, iso_utc

API = "https://ensemble-api.open-meteo.com/v1/ensemble"
CONFIG = RAIZ / "configs/evaluacion_final.jsonc"
LATITUD = -31.1115
LONGITUD = -71.6147


def cargar_config() -> dict:
    texto = pipe.quitar_comentarios_jsonc(CONFIG.read_text(encoding="utf-8"))
    return json.loads(texto)


def main() -> int:
    cfg = cargar_config()
    inicio = a_utc(cfg["inicio_obra_utc"])
    fin = inicio + timedelta(hours=float(cfg["tiempo_final_h"]))
    parametros = {
        "latitude": f"{LATITUD:.4f}",
        "longitude": f"{LONGITUD:.4f}",
        "hourly": "temperature_2m",
        "models": "ncep_gefs05",
        "timezone": "UTC",
        "forecast_days": 35,
        "temperature_unit": "celsius",
    }
    url = f"{API}?{urlencode(parametros)}"
    with urlopen(url, timeout=90) as respuesta:
        bruto = respuesta.read()
    datos = json.loads(bruto.decode("utf-8"))
    hourly = datos["hourly"]
    tiempos = tuple(a_utc(t) for t in hourly["time"])
    if not tiempos or tiempos[0] > inicio or tiempos[-1] < fin:
        raise ValueError(f"GEFS no cubre {iso_utc(inicio)} a {iso_utc(fin)}; "
                         f"cobertura: {iso_utc(tiempos[0])} a {iso_utc(tiempos[-1])}")
    indices = [i for i, t in enumerate(tiempos) if inicio <= t <= fin]
    tiempos_objetivo = tuple(tiempos[i] for i in indices)
    esperado = int(round((fin - inicio).total_seconds() / 3600.0)) + 1
    if (tiempos_objetivo[0] != inicio or tiempos_objetivo[-1] != fin
            or len(tiempos_objetivo) != esperado):
        raise ValueError("GEFS no entrega la grilla horaria exacta del horizonte solicitado")

    unidades = datos.get("hourly_units", {}).get("temperature_2m")
    if unidades not in ("°C", "degC"):
        raise ValueError(f"GEFS entregó una unidad inesperada para temperatura: {unidades!r}")
    claves = [("control", "temperature_2m")]
    claves.extend((f"member_{n:02d}", f"temperature_2m_member{n:02d}")
                  for n in range(1, 31))
    if any(clave not in hourly for _, clave in claves):
        raise ValueError("La respuesta GEFS no contiene el control y los 30 miembros esperados")
    for miembro, clave in claves:
        valores = [hourly[clave][i] for i in indices]
        if len(valores) != esperado or any(v is None or not np.isfinite(float(v)) for v in valores):
            raise ValueError(f"GEFS tiene temperaturas faltantes en {miembro} dentro del horizonte")

    recuperado = datetime.now(timezone.utc)
    marca = recuperado.strftime("%Y-%m-%dT%H%M%SZ")
    raiz_salida = RAIZ / "data/weather"
    raiz_salida.mkdir(parents=True, exist_ok=True)
    base = raiz_salida / f"{marca}_GEFS_ensemble"
    sufijo = 2
    while base.exists():
        base = raiz_salida / f"{marca}_GEFS_ensemble_{sufijo}"
        sufijo += 1
    base.mkdir()
    (base / "original.json").write_bytes(bruto)
    sha = hashlib.sha256(bruto).hexdigest()
    ids = []

    for miembro, clave in claves:
        carpeta = base / "members" / miembro
        carpeta.mkdir(parents=True)
        forecast = {
            "time_utc": [iso_utc(t) for t in tiempos_objetivo],
            "temperature_2m": [round(float(hourly[clave][i]), 4) for i in indices],
            "units": {"temperature_2m": "degC"},
        }
        (carpeta / "forecast.json").write_text(
            json.dumps(forecast, indent=2), encoding="utf-8")
        metadatos = {
            "source": "Open-Meteo Ensemble API; NOAA Global Ensemble Forecast System",
            "model": "ncep_gefs05",
            "member": miembro,
            "site": "Grupo 3: Mantos de Hornillos / Ovalle, Camarico – Ecos del Mar, Coquimbo",
            "latitude": LATITUD,
            "longitude": LONGITUD,
            "timezone": "UTC",
            "coordinates_source": "Tabla docente de sitios del curso; mismos datos usados en E11",
            "issued_at": None,
            "issued_at_note": "la API no expone la hora de emisión del ciclo GEFS",
            "retrieved_at": iso_utc(recuperado),
            "valid_from": iso_utc(inicio),
            "valid_to": iso_utc(fin),
            "variables": ["temperature_2m"],
            "units": {"temperature_2m": "degC"},
            "original_file": "../../original.json",
            "original_sha256": sha,
            "ensemble_manifest_file": "../../ensemble_metadata.json",
            "n_points": esperado,
            "series_kind": "pronostico",
            "synthetic": False,
            "ensemble_member_key": clave,
            "native_temporal_resolution": "3-hourly GEFS output; API time series hourly-interpolated",
        }
        (carpeta / "metadata.json").write_text(
            json.dumps(metadatos, indent=2, ensure_ascii=False), encoding="utf-8")
        ids.append({"id": miembro, "key": clave,
                    "snapshot": (carpeta / "metadata.json").parent.relative_to(RAIZ).as_posix()})

    resumen = {
        "source": "Open-Meteo Ensemble API",
        "model": "NOAA GEFS 0.5° (ncep_gefs05)",
        "retrieved_at": iso_utc(recuperado),
        "forecast_reference_start": iso_utc(tiempos[0]),
        "forecast_reference_end": iso_utc(tiempos[-1]),
        "target_start": iso_utc(inicio),
        "target_end": iso_utc(fin),
        "target_points_per_member": esperado,
        "native_grid_approx_km": 50,
        "native_time_step": "3 h",
        "api_time_step": "1 h (interpolated by Open-Meteo)",
        "member_count": len(claves),
        "member_weights": "equal, 1/31 each (provisional empirical ensemble weights)",
        "coordinates": {"latitude": LATITUD, "longitude": LONGITUD},
        "request_url": url,
        "original_file": "original.json",
        "original_sha256": sha,
        "members": ids,
    }
    (base / "ensemble_metadata.json").write_text(
        json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"Ensamble guardado: {base.relative_to(RAIZ)}")
    print(f"Cobertura usada: {iso_utc(inicio)} a {iso_utc(fin)} ({esperado} muestras por miembro)")
    print(f"Miembros recibidos y archivados: {len(claves)}")
    import preparar_clima_representativo_gefs as preparar
    preparar.main(base)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR DE DESCARGA/COBERTURA: {exc}", file=sys.stderr)
        raise SystemExit(2)
