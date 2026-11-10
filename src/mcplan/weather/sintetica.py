"""Serie de aire SINTÉTICA para demos y pruebas. No es un pronóstico ni una observación."""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from .serie import ORIGEN_SINTETICA, SerieTemperatura, a_utc

AVISO_SINTETICA = ("SERIE SINTÉTICA de demostración: media + ciclo diario sinusoidal + "
                   "tendencia lineal. No es un pronóstico ni un dato observado.")


def serie_sintetica(inicio_utc, duracion_h: float, *, media_c: float = 16.0,
                    amplitud_diaria_c: float = 6.0, hora_pico_utc: float = 18.0,
                    tendencia_c_por_dia: float = 0.0, paso_h: float = 1.0) -> SerieTemperatura:
    """T(t) = media + tendencia·t_días + A·cos(2π (hora_UTC − hora_pico)/24).

    Determinista (sin ruido). ``hora_pico_utc`` es la hora UTC del máximo diario.
    """
    if duracion_h <= 0 or paso_h <= 0:
        raise ValueError("duracion_h y paso_h deben ser positivos")
    inicio = a_utc(inicio_utc)
    n = int(np.floor(duracion_h / paso_h + 1e-9)) + 1
    horas_rel = np.arange(n) * paso_h
    tiempos = tuple(inicio + timedelta(hours=float(h)) for h in horas_rel)
    hora_dia = np.array([t.hour + t.minute / 60.0 for t in tiempos])
    temperatura = (media_c + tendencia_c_por_dia * horas_rel / 24.0
                   + amplitud_diaria_c * np.cos(2.0 * np.pi * (hora_dia - hora_pico_utc) / 24.0))
    return SerieTemperatura(tiempos, temperatura, tuple(ORIGEN_SINTETICA for _ in tiempos))
