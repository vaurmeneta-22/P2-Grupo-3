"""Modelos de fuente volumétrica dependientes de la edad de una vertida.

Incluye una curva adiabática exponencial y el pulso cosenoidal usado por el
pipeline IFC 3D. Sus parámetros son configurables y no son una calibración
del hormigón real del proyecto.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ModeloHidratacionExponencial:
    """Elevación adiabática asintótica y tiempo característico.

    ``T_ad(age) = T_colocación + delta_T * (1-exp(-age/tau))``.
    La fuente se deriva como ``Qdot = rho*cp*dT_ad/dt`` [W/m³].
    """

    incremento_adiabatico_c: float
    tiempo_caracteristico_s: float

    def __post_init__(self):
        if not math.isfinite(self.incremento_adiabatico_c) or self.incremento_adiabatico_c < 0:
            raise ValueError("incremento_adiabatico_c debe ser finito y no negativo")
        if not math.isfinite(self.tiempo_caracteristico_s) or self.tiempo_caracteristico_s <= 0:
            raise ValueError("tiempo_caracteristico_s debe ser positivo y finito")

    def temperatura_adiabatica_c(self, temperatura_colocacion_c: float,
                                 edad_s) -> np.ndarray | float:
        edad = np.maximum(np.asarray(edad_s, dtype=float), 0.0)
        temperatura = temperatura_colocacion_c + self.incremento_adiabatico_c * (
            1.0 - np.exp(-edad / self.tiempo_caracteristico_s))
        return float(temperatura) if temperatura.ndim == 0 else temperatura

    def generacion_w_m3(self, edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Devuelve Qdot [W/m³], nula antes de colocar y positiva después."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(rho <= 0) or np.any(cp <= 0):
            raise ValueError("rho y cp deben ser positivos")
        fuente = rho * cp * self.incremento_adiabatico_c / self.tiempo_caracteristico_s
        return np.where(edad >= 0.0, fuente * np.exp(-np.maximum(edad, 0.0)
                                                      / self.tiempo_caracteristico_s), 0.0)


@dataclass(frozen=True)
class ModeloPulsoCosenoidal:
    """Fuente volumétrica cosenoidal definida por la edad de cada vertida.

    ``Q(age) = Qmax/2 * (1-cos(2*pi*age/duracion))`` para
    ``0 <= age <= duracion`` y cero fuera de ese intervalo. La integral
    adiabática del pulso es ``Qmax*duracion/(2*rho*cp)`` [K].
    """

    qmax_w_m3: float
    duracion_s: float

    def __post_init__(self):
        if not math.isfinite(self.qmax_w_m3) or self.qmax_w_m3 < 0:
            raise ValueError("qmax_w_m3 debe ser finito y no negativo")
        if not math.isfinite(self.duracion_s) or self.duracion_s <= 0:
            raise ValueError("duracion_s debe ser positiva y finita")

    def generacion_w_m3(self, edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Devuelve Q(age) en W/m³; cada elemento puede tener edad distinta."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
            raise ValueError("rho debe ser finita y positiva")
        if np.any(~np.isfinite(cp)) or np.any(cp <= 0):
            raise ValueError("cp debe ser finito y positivo")
        activa = (edad >= 0.0) & (edad <= self.duracion_s)
        edad_segura = np.clip(edad, 0.0, self.duracion_s)
        q = 0.5 * self.qmax_w_m3 * (
            1.0 - np.cos(2.0 * math.pi * edad_segura / self.duracion_s)
        )
        return np.where(activa, q, 0.0)

    def temperatura_adiabatica_c(self, temperatura_colocacion_c: float,
                                  edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Temperatura acumulada ideal sin pérdidas para una vertida [°C]."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
            raise ValueError("rho debe ser finita y positiva")
        if np.any(~np.isfinite(cp)) or np.any(cp <= 0):
            raise ValueError("cp debe ser finito y positivo")
        edad_limitada = np.clip(edad, 0.0, self.duracion_s)
        energia_por_volumen = 0.5 * self.qmax_w_m3 * (
            edad_limitada - self.duracion_s / (2.0 * math.pi)
            * np.sin(2.0 * math.pi * edad_limitada / self.duracion_s)
        )
        temperatura = temperatura_colocacion_c + energia_por_volumen / (rho * cp)
        return float(temperatura) if temperatura.ndim == 0 else temperatura
