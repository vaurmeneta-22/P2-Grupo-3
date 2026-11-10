"""Crea el snapshot SINTÉTICO que usa la demo E11 mientras no haya coordenadas oficiales.

No es un pronóstico ni un dato observado. Es idempotente: si la carpeta ya existe no la toca.
Ejecutar con ▶️ o: python scripts/crear_snapshot_sintetico_e11.py
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from mcplan.weather.sintetica import AVISO_SINTETICA, serie_sintetica
from mcplan.weather.snapshots import guardar_snapshot

DESCARGA = datetime(2026, 11, 10, 16, 0, tzinfo=timezone.utc)      # fija el nombre de la carpeta
CARPETA = RAIZ / "data" / "weather" / "2026-11-10T160000Z"
INICIO = "2026-11-11T12:00:00Z"


def main() -> int:
    if CARPETA.exists():
        print(f"Ya existe (no se modifica): {CARPETA}")
        return 0
    serie = serie_sintetica(INICIO, 240.0, media_c=16.0, amplitud_diaria_c=7.0,
                            hora_pico_utc=18.0, tendencia_c_por_dia=0.3)
    bruto = json.dumps({"aviso": AVISO_SINTETICA, "inicio_utc": INICIO, "duracion_h": 240.0,
                        "parametros": {"media_c": 16.0, "amplitud_diaria_c": 7.0, "hora_pico_utc": 18.0,
                                       "tendencia_c_por_dia": 0.3}}, indent=2,
                       ensure_ascii=False).encode("utf-8")
    carpeta = guardar_snapshot(
        RAIZ / "data" / "weather", bruto=bruto, nombre_original="original.json", serie=serie,
        fuente="serie sintética (scripts/crear_snapshot_sintetico_e11.py)", modelo="sintético",
        latitud=None, longitud=None, descargado_utc=DESCARGA, emitido_utc=None,
        sitio="Grupo 3: Ovalle / Ecos del Mar, Coquimbo (coordenadas pendientes)", tipo="sintetica")
    print(f"Snapshot sintético creado: {carpeta}\n{AVISO_SINTETICA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
