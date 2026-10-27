"""Soluciones exactas y verificaciones reutilizables (tests y scripts/verificar_E05.py).

Todas las temperaturas en °C, longitudes en m, tiempos en s, k en W/(m K),
rho en kg/m³, cp en J/(kg K), q0 en W/m³ y energías en J. Las soluciones exactas
se evalúan en los CENTROS de celda: x_i = (i + 1/2)·Lx/nx.
"""
from dataclasses import dataclass
import math

import numpy as np

from mcplan.thermal.diffusion_nd import DifusionCartesiana, aislado, dirichlet

CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")


def centros_de_celda(forma, longitudes_m):
    """Coordenadas de centros de celda por eje (arreglos 1D)."""
    return [(np.arange(n) + 0.5) * (L / n) for n, L in zip(forma, longitudes_m)]


def lambda_modo(longitudes_m) -> float:
    """λ = Σ (π/L_eje)² [1/m²]."""
    return sum((math.pi / L) ** 2 for L in longitudes_m)


def solucion_sinusoidal(forma, longitudes_m, T0_c, amplitud_c, alfa_m2_s, tiempo_s):
    """T = T0 + A·Π sin(π x_e/L_e)·exp(−α λ t) en los centros de celda."""
    campo = np.ones(forma)
    for m, L in zip(np.meshgrid(*centros_de_celda(forma, longitudes_m), indexing="ij"),
                    longitudes_m):
        campo = campo * np.sin(math.pi * m / L)
    return T0_c + amplitud_c * campo * math.exp(-alfa_m2_s * lambda_modo(longitudes_m) * tiempo_s)


@dataclass(frozen=True)
class ResultadoVerificacion:
    nombre: str
    error_max_c: float
    tolerancia_c: float
    detalle: dict

    @property
    def aprobado(self) -> bool:
        return bool(self.error_max_c <= self.tolerancia_c)


def verificar_sinusoidal(*, forma, longitudes_m, T0_c, amplitud_c, k_w_m_k, rho_kg_m3,
                         cp_j_kg_k, tiempo_final_s, fraccion_paso=0.9, paso_tiempo_max_s=None,
                         constante_tolerancia=1.0):
    """Dirichlet T0 en todas las caras, Qdot = 0, condición inicial exacta en t = 0.

    Con paso_tiempo_max_s=None usa fraccion_paso·dt_estable. El criterio usa el dt
    EFECTIVO de la corrida (ver ``estimar_error_discretizacion``). amplitud_c = 0 es
    válido (equilibrio): el criterio queda en un piso absoluto de redondeo.
    """
    ndim = len(forma)
    dx = tuple(L / n for L, n in zip(longitudes_m, forma))
    s = DifusionCartesiana(forma=forma, dx_m=dx, conductividad_w_m_k=k_w_m_k,
                           densidad_kg_m3=rho_kg_m3, calor_especifico_j_kg_k=cp_j_kg_k,
                           fronteras={c: dirichlet(T0_c) for c in CARAS[:2 * ndim]})
    alfa = k_w_m_k / (rho_kg_m3 * cp_j_kg_k)
    T_ini = solucion_sinusoidal(forma, longitudes_m, T0_c, amplitud_c, alfa, 0.0)
    dt = paso_tiempo_max_s or fraccion_paso * s.paso_maximo_estable_s()
    r = s.simular(temperatura_inicial_c=T_ini, tiempo_final_s=tiempo_final_s, paso_tiempo_max_s=dt)
    exacta = solucion_sinusoidal(forma, longitudes_m, T0_c, amplitud_c, alfa, tiempo_final_s)
    error = float(np.max(np.abs(r.temperatura_c - exacta)))
    est_esp, est_tmp = estimar_error_discretizacion(longitudes_m, dx, amplitud_c, alfa,
                                                    tiempo_final_s, r.paso_tiempo_s)
    piso = 1e-12 * max(1.0, abs(T0_c) + abs(amplitud_c))       # redondeo de maquina [°C]
    tol = constante_tolerancia * (est_esp + est_tmp) + piso
    return ResultadoVerificacion(
        f"sinusoidal {ndim}D", error, tol,
        {"error_normalizado": error / abs(amplitud_c) if amplitud_c != 0 else None,
         "error_sobre_tolerancia": error / tol,
         "estimacion_espacial_c": est_esp, "estimacion_temporal_c": est_tmp,
         "dt_s": r.paso_tiempo_s, "dt_limite_s": r.paso_maximo_estable_s, "pasos": r.numero_pasos,
         "factor_decaimiento": math.exp(-alfa * lambda_modo(longitudes_m) * tiempo_final_s),
         "T_max_numerica_c": float(r.temperatura_c.max()), "T_max_exacta_c": float(exacta.max()),
         "error_balance_j": r.error_balance_j, "alfa_m2_s": alfa})


def estimar_error_discretizacion(longitudes_m, dx_m, amplitud_c, alfa, tiempo_s, dt_s):
    """ESTIMACIÓN asintótica (no cota garantizada) del error máximo [°C]: (espacial, temporal).

    Sea μ = α λ la tasa exacta de decaimiento del modo, y e_rel el error relativo
    de la amplitud en el instante t. A primer orden:
    * Espacial: cada eje tiene autovalor discreto (4/dx²) sin²(π dx/2L) =
      (π/L)²[1 − (π dx/L)²/12 + ...], así que la tasa numérica difiere en
      α Δλ, Δλ = Σ_e (π/L_e)⁴ dx_e²/12  ⇒  e_esp = |A| e^{−μt} α t Δλ.
    * Temporal (Euler adelantado, dt EFECTIVO): (1 − μ dt)^n vs e^{−μ t} difiere en
      e^{−μt} μ² dt t/2  ⇒  e_tmp = |A| e^{−μt} μ² dt t/2.
    Los dos tienen signo opuesto (el espacial frena el decaimiento, Euler lo acelera),
    por lo que pueden compensarse; el criterio usa la SUMA DE MAGNITUDES para no
    depender de esa compensación. No incluye términos de orden superior ni el efecto
    de la condición de borde de media celda: por eso es una estimación a validar
    contra mallas y pasos distintos (tests), no una cota rigurosa.
    """
    mu = alfa * lambda_modo(longitudes_m)
    d_lam = sum((math.pi / L) ** 4 * d * d / 12.0 for d, L in zip(dx_m, longitudes_m))
    base = abs(amplitud_c) * math.exp(-mu * tiempo_s)
    return base * alfa * tiempo_s * d_lam, base * mu * mu * dt_s * tiempo_s / 2.0


def verificar_fuente_uniforme(*, forma, dx_m, T_inicial_c, q0_w_m3, k_w_m_k, rho_kg_m3,
                              cp_j_kg_k, tiempo_final_s, tolerancia_relativa=1e-10):
    """Dominio adiabático con Qdot uniforme: T = T_ini + q0 t/(rho cp) y
    E_final − E_inicial = q0·V_total·t (V_total = volumen de todas las celdas).

    q0 = 0 es válido (equilibrio). Criterios (redondeo de máquina, no truncamiento,
    pues la fuente uniforme es exacta en Euler): |ΔT| ≤ tol_rel·max(|T_esp|, 1) y
    |ΔE − ΔE_esp| ≤ tol_rel·max(|ΔE_esp|, ρ cp V_total·max(|T_ini|, 1)); el segundo
    denominador evita dividir entre ΔE_esp = 0 usando la energía interna como escala.
    """
    ndim = len(forma)
    s = DifusionCartesiana(forma=forma, dx_m=dx_m, conductividad_w_m_k=k_w_m_k,
                           densidad_kg_m3=rho_kg_m3, calor_especifico_j_kg_k=cp_j_kg_k,
                           fuente_w_m3=q0_w_m3,
                           fronteras={c: aislado() for c in CARAS[:2 * ndim]})
    r = s.simular(temperatura_inicial_c=T_inicial_c, tiempo_final_s=tiempo_final_s,
                  paso_tiempo_max_s=s.paso_maximo_estable_s())
    T_esp = T_inicial_c + q0_w_m3 * tiempo_final_s / (rho_kg_m3 * cp_j_kg_k)
    error_T = float(np.max(np.abs(r.temperatura_c - T_esp)))
    V_total = s.volumen_m3 * float(np.prod(forma))
    dE_esp = q0_w_m3 * V_total * tiempo_final_s
    dE = r.energia_final_j - r.energia_inicial_j
    escala_j = rho_kg_m3 * cp_j_kg_k * V_total * max(abs(T_inicial_c), 1.0)
    denominador = max(abs(dE_esp), escala_j)
    tol_T = tolerancia_relativa * max(abs(T_esp), 1.0)
    return ResultadoVerificacion(
        f"fuente uniforme {ndim}D", error_T, tol_T,
        {"T_esperada_c": T_esp, "dE_j": dE, "dE_esperado_j": dE_esp,
         "error_energia_j": abs(dE - dE_esp), "denominador_energia_j": denominador,
         "error_energia_relativo": abs(dE - dE_esp) / denominador,
         "tol_energia_relativa": tolerancia_relativa,
         "V_total_m3": V_total, "dt_s": r.paso_tiempo_s, "pasos": r.numero_pasos})
