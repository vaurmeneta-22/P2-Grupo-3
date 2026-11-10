"""Descarga un snapshot REAL de Open-Meteo para el sitio del Grupo 3.

Exige las coordenadas oficiales de sitios.yaml: en ``configs/e11_clima_real.jsonc``
(bloque ``sitio``) o con --lat/--lon. Sin coordenadas termina con un mensaje claro y no
toca la red. Cada descarga crea una carpeta nueva en data/weather/ (no sobrescribe).

  Pronóstico (hasta 16 días):    python scripts/descargar_snapshot_e11.py --modo pronostico
  Serie histórica (reanálisis):  python scripts/descargar_snapshot_e11.py --modo historico --inicio 2025-11-11 --fin 2025-11-25

Después, poner la carpeta impresa en ``clima.snapshot.carpeta`` de la configuración E11.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from mcplan.weather.openmeteo import descargar_historica_y_guardar, descargar_y_guardar
from mcplan.weather.sitio import SitioSinCoordenadas, sitio_desde_config

CONFIG = RAIZ / "configs" / "e11_clima_real.jsonc"


def quitar_comentarios(texto: str) -> str:
    salida, en_cadena, escape, i = [], False, False, 0
    while i < len(texto):
        c = texto[i]
        if en_cadena:
            salida.append(c)
            escape = (c == "\\" and not escape)
            if c == '"' and not escape:
                en_cadena = False
        elif c == '"':
            en_cadena = True
            salida.append(c)
        elif texto.startswith("//", i):
            while i < len(texto) and texto[i] != "\n":
                i += 1
            continue
        else:
            salida.append(c)
        i += 1
    return "".join(salida)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--modo", choices=("pronostico", "historico"), default="pronostico")
    parser.add_argument("--inicio", help="Solo histórico: YYYY-MM-DD")
    parser.add_argument("--fin", help="Solo histórico: YYYY-MM-DD")
    parser.add_argument("--lat", type=float, help="Sustituye la latitud del config")
    parser.add_argument("--lon", type=float, help="Sustituye la longitud del config")
    parser.add_argument("--modelo", help="Producto/modelo de Open-Meteo (opcional)")
    args = parser.parse_args()
    config = json.loads(quitar_comentarios(CONFIG.read_text(encoding="utf-8")))
    cfg_sitio = dict(config.get("sitio") or {})
    if args.lat is not None:
        cfg_sitio["latitud"] = args.lat
    if args.lon is not None:
        cfg_sitio["longitud"] = args.lon
    sitio = sitio_desde_config(cfg_sitio)
    raiz = RAIZ / "data" / "weather"
    try:
        if args.modo == "pronostico":
            carpeta = descargar_y_guardar(sitio, raiz, modelo=args.modelo)
        else:
            if not (args.inicio and args.fin):
                parser.error("--modo historico requiere --inicio y --fin")
            carpeta = descargar_historica_y_guardar(sitio, raiz, args.inicio, args.fin, modelo=args.modelo)
    except SitioSinCoordenadas as error:
        print(f"DESCARGA BLOQUEADA (no se consultó la red): {error}")
        return 2
    print(f"Snapshot guardado: {carpeta}\nPonlo en clima.snapshot.carpeta de configs/e11_clima_real.jsonc "
          f"como \"data/weather/{carpeta.name}\".")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
