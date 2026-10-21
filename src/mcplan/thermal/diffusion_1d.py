"""Solver didáctico 1D de difusión con bordes Dirichlet, Neumann y Robin.

Implementa FTCS explícito. En Neumann se prescribe el gradiente dT/dx,
en K/m; gradiente cero representa un extremo aislado. Robin modela
convección hacia un ambiente a temperatura prescrita. También permite una
generación volumétrica uniforme de calor, en W/m³.
"""
from dataclasses import dataclass
from collections.abc import Callable
import math

import numpy as np


@dataclass(frozen=True)
class ResultadoDifusion:
    """Resultado del último instante de una simulación 1D."""

    posicion_m: np.ndarray
    temperatura_c: np.ndarray
    paso_tiempo_s: float
    numero_pasos: int
    numero_fourier: float


def calcular_difusividad_termica(conductividad_w_m_k: float,
                                 densidad_kg_m3: float,
                                 calor_especifico_j_kg_k: float) -> float:
    """Calcula alfa = k/(rho*cp), en m²/s, con propiedades SI."""
    propiedades = (conductividad_w_m_k, densidad_kg_m3,
                   calor_especifico_j_kg_k)
    if not all(math.isfinite(valor) and valor > 0 for valor in propiedades):
        raise ValueError("k, rho y cp deben ser valores finitos y positivos")
    return conductividad_w_m_k / (densidad_kg_m3 * calor_especifico_j_kg_k)


def simular_difusion_1d(
    *,
    temperatura_inicial_c,
    longitud_m: float,
    difusividad_m2_s: float,
    condicion_izquierda: str,
    valor_izquierdo: float | None,
    condicion_derecha: str,
    valor_derecho: float | None,
    conductividad_w_m_k: float | None = None,
    coeficiente_conveccion_izquierdo_w_m2_k: float = 0.0,
    coeficiente_conveccion_derecho_w_m2_k: float = 0.0,
    temperatura_ambiente_izquierda_c: float | Callable[[float], float] | None = None,
    temperatura_ambiente_derecha_c: float | Callable[[float], float] | None = None,
    generacion_calor_w_m3: float = 0.0,
    capacidad_volumetrica_j_m3_k: float | None = None,
    tiempo_final_s: float,
    paso_tiempo_max_s: float,
) -> ResultadoDifusion:
    """Resuelve difusión 1D con FTCS y condiciones independientes por extremo.

    ``condicion_*`` admite ``"dirichlet"`` (valor en °C), ``"neumann"``
    (gradiente dT/dx en K/m) o ``"robin"``. Para Robin, la temperatura
    ambiente puede ser un valor constante o una función del tiempo en segundos;
    también se requieren conductividad ``k`` y coeficientes ``h`` en W/(m² K).
    ``generacion_calor_w_m3`` agrega una fuente uniforme Q₀. Si Q₀ no es cero,
    también se requiere ``capacidad_volumetrica_j_m3_k`` (rho*cp). Se verifica
    la estabilidad, incluyendo el efecto de Robin en los bordes.
    """
    condicion_izquierda = condicion_izquierda.lower()
    condicion_derecha = condicion_derecha.lower()
    condiciones_validas = {"dirichlet", "neumann", "robin"}
    if condicion_izquierda not in condiciones_validas or condicion_derecha not in condiciones_validas:
        raise ValueError("Las condiciones deben ser 'dirichlet', 'neumann' o 'robin'")

    parametros = (longitud_m, difusividad_m2_s, tiempo_final_s, paso_tiempo_max_s,
                  coeficiente_conveccion_izquierdo_w_m2_k,
                  coeficiente_conveccion_derecho_w_m2_k,
                  generacion_calor_w_m3)
    if not all(math.isfinite(valor) for valor in parametros):
        raise ValueError("Los parámetros deben ser finitos")
    if (longitud_m <= 0 or difusividad_m2_s <= 0 or tiempo_final_s <= 0
            or paso_tiempo_max_s <= 0):
        raise ValueError("L, alfa, tiempo final y paso deben ser positivos")
    if (coeficiente_conveccion_izquierdo_w_m2_k < 0
            or coeficiente_conveccion_derecho_w_m2_k < 0):
        raise ValueError("Los coeficientes de convección h no pueden ser negativos")
    if generacion_calor_w_m3 != 0.0:
        if (capacidad_volumetrica_j_m3_k is None
                or not math.isfinite(capacidad_volumetrica_j_m3_k)
                or capacidad_volumetrica_j_m3_k <= 0):
            raise ValueError(
                "Q0 distinto de cero requiere rho*cp finito y positivo "
                "en J/(m³ K)"
            )

    for condicion, valor, nombre in (
        (condicion_izquierda, valor_izquierdo, "izquierdo"),
        (condicion_derecha, valor_derecho, "derecho"),
    ):
        if condicion in {"dirichlet", "neumann"}:
            if valor is None or not math.isfinite(valor):
                raise ValueError(
                    f"El valor del borde {nombre} debe ser numérico y finito"
                )

    usa_robin = "robin" in (condicion_izquierda, condicion_derecha)
    if usa_robin:
        if conductividad_w_m_k is None or not math.isfinite(conductividad_w_m_k):
            raise ValueError("Robin requiere una conductividad k finita y positiva")
        if conductividad_w_m_k <= 0:
            raise ValueError("Robin requiere una conductividad k finita y positiva")
        if condicion_izquierda == "robin":
            _evaluar_temperatura_ambiente(
                temperatura_ambiente_izquierda_c, 0.0, "izquierda"
            )
        if condicion_derecha == "robin":
            _evaluar_temperatura_ambiente(
                temperatura_ambiente_derecha_c, 0.0, "derecha"
            )
    try:
        temperatura = np.asarray(temperatura_inicial_c, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("La temperatura inicial debe ser un arreglo numérico 1D") from exc
    if temperatura.ndim != 1 or temperatura.size < 3:
        raise ValueError("Se requieren al menos tres temperaturas en un arreglo 1D")
    if not np.all(np.isfinite(temperatura)):
        raise ValueError("La temperatura inicial debe contener valores finitos")

    numero_nodos = temperatura.size
    posicion = np.linspace(0.0, longitud_m, numero_nodos)
    paso_espacial_m = posicion[1] - posicion[0]

    biot_malla_izquierdo = (
        coeficiente_conveccion_izquierdo_w_m2_k * paso_espacial_m
        / conductividad_w_m_k
        if condicion_izquierda == "robin" else 0.0
    )
    biot_malla_derecho = (
        coeficiente_conveccion_derecho_w_m2_k * paso_espacial_m
        / conductividad_w_m_k
        if condicion_derecha == "robin" else 0.0
    )
    fourier_maximo_estable = min(
        0.5,
        0.5 / (1.0 + biot_malla_izquierdo),
        0.5 / (1.0 + biot_malla_derecho),
    )
    numero_pasos = math.ceil(tiempo_final_s / paso_tiempo_max_s)
    paso_tiempo_s = tiempo_final_s / numero_pasos
    numero_fourier = difusividad_m2_s * paso_tiempo_s / paso_espacial_m**2
    incremento_fuente_c = (
        paso_tiempo_s * generacion_calor_w_m3 / capacidad_volumetrica_j_m3_k
        if generacion_calor_w_m3 != 0.0 else 0.0
    )
    if numero_fourier > fourier_maximo_estable + 1e-14:
        raise ValueError(
            "Paso inestable para FTCS 1D: se requiere "
            f"alfa*dt/dx² <= {fourier_maximo_estable:.6g}; "
            f"se obtuvo {numero_fourier:.6g}"
        )

    temperatura = temperatura.copy()
    if condicion_izquierda == "dirichlet":
        temperatura[0] = valor_izquierdo
    if condicion_derecha == "dirichlet":
        temperatura[-1] = valor_derecho

    for numero_paso in range(numero_pasos):
        tiempo_actual_s = numero_paso * paso_tiempo_s
        temperatura_anterior = temperatura
        temperatura_siguiente = temperatura_anterior.copy()
        temperatura_siguiente[1:-1] = temperatura_anterior[1:-1] + numero_fourier * (
            temperatura_anterior[2:] - 2.0 * temperatura_anterior[1:-1]
            + temperatura_anterior[:-2]
        ) + incremento_fuente_c
        if condicion_izquierda == "dirichlet":
            temperatura_siguiente[0] = valor_izquierdo
        elif condicion_izquierda == "neumann":
            # Punto fantasma: dT/dx = (T_1 - T_-1)/(2 dx).
            temperatura_siguiente[0] = temperatura_anterior[0] + 2.0 * numero_fourier * (
                temperatura_anterior[1] - temperatura_anterior[0]
                - paso_espacial_m * valor_izquierdo
            ) + incremento_fuente_c
        else:
            # Balance en el medio volumen del borde: conducción + convección.
            temperatura_aire_izquierda_c = _evaluar_temperatura_ambiente(
                temperatura_ambiente_izquierda_c, tiempo_actual_s, "izquierda"
            )
            temperatura_siguiente[0] = temperatura_anterior[0] + 2.0 * numero_fourier * (
                temperatura_anterior[1] - temperatura_anterior[0]
                + biot_malla_izquierdo * (
                    temperatura_aire_izquierda_c - temperatura_anterior[0]
                )
            ) + incremento_fuente_c
        if condicion_derecha == "dirichlet":
            temperatura_siguiente[-1] = valor_derecho
        elif condicion_derecha == "neumann":
            # Punto fantasma: dT/dx = (T_N - T_(N-2))/(2 dx).
            temperatura_siguiente[-1] = temperatura_anterior[-1] + 2.0 * numero_fourier * (
                temperatura_anterior[-2] - temperatura_anterior[-1]
                + paso_espacial_m * valor_derecho
            ) + incremento_fuente_c
        else:
            temperatura_aire_derecha_c = _evaluar_temperatura_ambiente(
                temperatura_ambiente_derecha_c, tiempo_actual_s, "derecha"
            )
            temperatura_siguiente[-1] = temperatura_anterior[-1] + 2.0 * numero_fourier * (
                temperatura_anterior[-2] - temperatura_anterior[-1]
                + biot_malla_derecho * (
                    temperatura_aire_derecha_c - temperatura_anterior[-1]
                )
            ) + incremento_fuente_c
        temperatura = temperatura_siguiente

    return ResultadoDifusion(
        posicion_m=posicion,
        temperatura_c=temperatura,
        paso_tiempo_s=paso_tiempo_s,
        numero_pasos=numero_pasos,
        numero_fourier=numero_fourier,
    )


def _evaluar_temperatura_ambiente(
    temperatura_ambiente_c: float | Callable[[float], float] | None,
    tiempo_s: float,
    lado: str,
) -> float:
    """Evalúa una temperatura ambiente constante o dependiente del tiempo."""
    if temperatura_ambiente_c is None:
        raise ValueError(f"Robin requiere la temperatura ambiente del borde {lado}")
    try:
        valor = (
            temperatura_ambiente_c(tiempo_s)
            if callable(temperatura_ambiente_c)
            else temperatura_ambiente_c
        )
        valor = float(valor)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"La temperatura ambiente del borde {lado} debe ser numérica"
        ) from exc
    if not math.isfinite(valor):
        raise ValueError(
            f"La temperatura ambiente del borde {lado} debe ser finita"
        )
    return valor


def simular_dirichlet(
    *,
    temperatura_inicial_c,
    longitud_m: float,
    difusividad_m2_s: float,
    temperatura_izquierda_c: float,
    temperatura_derecha_c: float,
    tiempo_final_s: float,
    paso_tiempo_max_s: float,
) -> ResultadoDifusion:
    """Compatibilidad: resuelve el caso Dirichlet en ambos extremos."""
    return simular_difusion_1d(
        temperatura_inicial_c=temperatura_inicial_c,
        longitud_m=longitud_m,
        difusividad_m2_s=difusividad_m2_s,
        condicion_izquierda="dirichlet",
        valor_izquierdo=temperatura_izquierda_c,
        condicion_derecha="dirichlet",
        valor_derecho=temperatura_derecha_c,
        tiempo_final_s=tiempo_final_s,
        paso_tiempo_max_s=paso_tiempo_max_s,
    )
