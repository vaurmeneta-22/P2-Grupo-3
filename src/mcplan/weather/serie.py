"""Serie temporal de temperatura del aire con la procedencia de cada punto."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Callable

import numpy as np

ORIGEN_PRONOSTICO = "pronostico"
ORIGEN_SINTETICA = "sintetica"


def a_utc(valor) -> datetime:
    """Convierte texto ISO o ``datetime`` a ``datetime`` consciente de zona en UTC.

    Un valor sin zona se interpreta como UTC (los snapshots guardan siempre UTC).
    """
    if isinstance(valor, str):
        valor = datetime.fromisoformat(valor.replace("Z", "+00:00"))
    if valor.tzinfo is None:
        valor = valor.replace(tzinfo=timezone.utc)
    return valor.astimezone(timezone.utc)


def iso_utc(valor: datetime) -> str:
    return a_utc(valor).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass(frozen=True)
class SerieTemperatura:
    """Temperatura del aire [°C] en instantes UTC crecientes.

    ``origen[i]`` distingue el pronóstico original (``"pronostico"``), la serie
    sintética (``"sintetica"``) y cualquier extensión (``"extension:<método>"``).
    """

    tiempos_utc: tuple[datetime, ...]
    temperatura_c: np.ndarray
    origen: tuple[str, ...]

    def __post_init__(self):
        n = len(self.tiempos_utc)
        if n < 2:
            raise ValueError("la serie necesita al menos dos puntos")
        if len(self.temperatura_c) != n or len(self.origen) != n:
            raise ValueError("tiempos, temperaturas y orígenes deben tener la misma longitud")
        if not np.all(np.isfinite(self.temperatura_c)):
            raise ValueError("la temperatura debe ser finita (sin NaN)")
        segundos = np.array([t.timestamp() for t in self.tiempos_utc])
        if np.any(np.diff(segundos) <= 0):
            raise ValueError("los tiempos deben ser estrictamente crecientes")

    @property
    def intervalo(self) -> tuple[datetime, datetime]:
        return self.tiempos_utc[0], self.tiempos_utc[-1]

    def segundos_desde(self, t0_utc: datetime) -> np.ndarray:
        t0 = a_utc(t0_utc)
        return np.array([(t - t0).total_seconds() for t in self.tiempos_utc])

    def a_funcion(self, t0_utc: datetime, *, permitir_fuera: bool = False) -> Callable[[float], float]:
        """Función ``T(t_s)`` con ``t_s`` medido desde ``t0_utc`` (interpolación lineal).

        Fuera de la cobertura lanza ``ValueError`` (nunca inventa datos), salvo con
        ``permitir_fuera=True``, que mantiene el valor del extremo.
        """
        x = self.segundos_desde(t0_utc)
        y = np.asarray(self.temperatura_c, dtype=float)

        def temperatura(tiempo_s: float) -> float:
            if not permitir_fuera and not (x[0] - 1e-6 <= tiempo_s <= x[-1] + 1e-6):
                raise ValueError(
                    f"t={tiempo_s / 3600:.2f} h está fuera de la serie "
                    f"[{x[0] / 3600:.2f}, {x[-1] / 3600:.2f}] h desde el inicio de obra")
            return float(np.interp(tiempo_s, x, y))

        return temperatura

    def fraccion_por_origen(self) -> dict[str, int]:
        conteo: dict[str, int] = {}
        for origen in self.origen:
            conteo[origen] = conteo.get(origen, 0) + 1
        return conteo


def serie_desde_listas(tiempos, temperaturas, origen: str) -> SerieTemperatura:
    return SerieTemperatura(tuple(a_utc(t) for t in tiempos),
                            np.asarray(temperaturas, dtype=float),
                            tuple(origen for _ in tiempos))


def horas(n: float) -> timedelta:
    return timedelta(hours=n)


TIPOS_SERIE = ("pronostico", "reanalisis", "observacion", "sintetica")


def normalizar_serie(tiempos, valores, *, unidad: str = "degC", origen: str = ORIGEN_PRONOSTICO,
                     faltantes: str = "error") -> SerieTemperatura:
    """Valida y normaliza una serie cruda a °C y UTC.

    Rechaza: unidades desconocidas (se aceptan degC, °C, K), tiempos duplicados o
    desordenados (no se reordena en silencio) y largos distintos. Los datos faltantes
    (None o NaN) generan error, o se descartan con ``faltantes="descartar"``
    (la cobertura real se informa luego con ``snapshots.cobertura``).
    """
    if len(tiempos) != len(valores):
        raise ValueError("tiempos y valores tienen distinto largo")
    if faltantes not in ("error", "descartar"):
        raise ValueError("faltantes debe ser 'error' o 'descartar'")
    unidad_norm = unidad.strip().lower()
    if unidad_norm in ("degc", "°c", "c", "celsius"):
        corregir = 0.0
    elif unidad_norm in ("k", "kelvin"):
        corregir = -273.15
    else:
        raise ValueError(f"unidad no soportada: {unidad!r} (use degC o K)")
    pares = []
    for t, v in zip(tiempos, valores):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            if faltantes == "error":
                raise ValueError(f"dato faltante en {t}")
            continue
        pares.append((a_utc(t), float(v) + corregir))
    marcas = [p[0] for p in pares]
    if len(set(marcas)) != len(marcas):
        raise ValueError("la serie tiene tiempos duplicados")
    if any(b <= a for a, b in zip(marcas, marcas[1:])):
        raise ValueError("los tiempos de la serie no están en orden creciente")
    return SerieTemperatura(tuple(marcas), np.array([p[1] for p in pares]),
                            tuple(origen for _ in pares))
