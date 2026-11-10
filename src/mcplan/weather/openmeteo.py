"""Primer cliente meteorologico: Open-Meteo (temperature_2m horaria).

* ``construir_url`` y ``parsear_respuesta`` no usan la red y se prueban con una
  respuesta guardada.
* ``descargar_y_guardar`` exige coordenadas verificadas del sitio ANTES de cualquier
  consulta (``SitioSinCoordenadas``): la falta de latitud/longitud solo bloquea la
  consulta real, no el resto del pipeline.
* El endpoint de pronostico no informa la hora de emision del modelo; queda
  ``issued_at = null`` salvo que se entregue ``emitido_utc``.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from typing import Callable
from urllib.parse import urlencode
from urllib.request import urlopen

from .serie import ORIGEN_PRONOSTICO, SerieTemperatura, normalizar_serie
from .sitio import Sitio
from .snapshots import guardar_snapshot

URL_PRONOSTICO = "https://api.open-meteo.com/v1/forecast"
URL_HISTORICO = "https://archive-api.open-meteo.com/v1/archive"   # reanalisis (ERA5 y similares)


def construir_url(latitud: float, longitud: float, *, dias: int = 16, modelo: str | None = None,
                  base: str = URL_PRONOSTICO) -> str:
    """URL de pronostico horario en UTC (``timezone=UTC`` evita ambiguedades de hora local)."""
    parametros = {"latitude": f"{latitud:.4f}", "longitude": f"{longitud:.4f}",
                  "hourly": "temperature_2m", "timezone": "UTC", "forecast_days": int(dias),
                  "temperature_unit": "celsius"}
    if modelo:
        parametros["models"] = modelo
    return f"{base}?{urlencode(parametros)}"


def construir_url_historica(latitud: float, longitud: float, inicio: str, fin: str, *,
                            modelo: str | None = None, base: str = URL_HISTORICO) -> str:
    """URL de la serie historica horaria (``inicio``/``fin`` como YYYY-MM-DD, en UTC)."""
    parametros = {"latitude": f"{latitud:.4f}", "longitude": f"{longitud:.4f}",
                  "start_date": inicio, "end_date": fin, "hourly": "temperature_2m",
                  "timezone": "UTC", "temperature_unit": "celsius"}
    if modelo:
        parametros["models"] = modelo
    return f"{base}?{urlencode(parametros)}"


def parsear_respuesta(bruto: bytes) -> SerieTemperatura:
    """Normaliza la respuesta JSON a una serie en UTC. Rechaza unidades distintas de grados C."""
    datos = json.loads(bruto.decode("utf-8"))
    if datos.get("hourly_units", {}).get("temperature_2m") not in ("°C", "degC"):
        raise ValueError("Open-Meteo: la temperatura no esta en grados Celsius")
    tiempos = datos["hourly"]["time"]
    valores = datos["hourly"]["temperature_2m"]
    if len(tiempos) != len(valores):
        raise ValueError("Open-Meteo: time y temperature_2m tienen distinto largo")
    if int(datos.get("utc_offset_seconds", 0)) != 0:
        raise ValueError("Open-Meteo: se esperaba timezone=UTC (utc_offset_seconds=0)")
    return normalizar_serie(tiempos, valores, unidad="degC", origen=ORIGEN_PRONOSTICO,
                            faltantes="descartar")


def _obtener(url: str) -> bytes:  # pragma: no cover (usa la red)
    with urlopen(url, timeout=30) as respuesta:
        return respuesta.read()


def descargar_y_guardar(sitio: Sitio, raiz, *, dias: int = 16, modelo: str | None = None,
                        emitido_utc: datetime | None = None,
                        obtener: Callable[[str], bytes] = _obtener,
                        ahora_utc: datetime | None = None):
    """Descarga el pronostico del sitio y crea un snapshot nuevo. Devuelve su carpeta."""
    latitud, longitud = sitio.exigir_coordenadas()          # antes de tocar la red
    url = construir_url(latitud, longitud, dias=dias, modelo=modelo)
    bruto = obtener(url)
    serie = parsear_respuesta(bruto)
    return guardar_snapshot(
        raiz, bruto=bruto, nombre_original="original.json", serie=serie,
        fuente="Open-Meteo forecast API", modelo=modelo or "best_match (por defecto de la API)",
        latitud=latitud, longitud=longitud, zona_horaria="UTC",
        descargado_utc=ahora_utc or datetime.now(timezone.utc), emitido_utc=emitido_utc,
        sitio=f"Grupo {sitio.grupo}: {sitio.nombre}, {sitio.region}",
        procedencia_coordenadas=sitio.fuente_coordenadas,
        extra={"url": url})


def descargar_historica_y_guardar(sitio: Sitio, raiz, inicio: str, fin: str, *,
                                  modelo: str | None = None,
                                  obtener: Callable[[str], bytes] = _obtener,
                                  ahora_utc: datetime | None = None):
    """Descarga una serie historica (tipo ``reanalisis``) y crea un snapshot nuevo."""
    latitud, longitud = sitio.exigir_coordenadas()          # antes de tocar la red
    url = construir_url_historica(latitud, longitud, inicio, fin, modelo=modelo)
    bruto = obtener(url)
    serie = parsear_respuesta(bruto)
    return guardar_snapshot(
        raiz, bruto=bruto, nombre_original="original.json", serie=serie,
        fuente="Open-Meteo historical weather API", modelo=modelo or "por defecto de la API",
        latitud=latitud, longitud=longitud, zona_horaria="UTC",
        descargado_utc=ahora_utc or datetime.now(timezone.utc), emitido_utc=None,
        sitio=f"Grupo {sitio.grupo}: {sitio.nombre}, {sitio.region}",
        procedencia_coordenadas=sitio.fuente_coordenadas, tipo="reanalisis",
        extra={"url": url, "inicio": inicio, "fin": fin})
