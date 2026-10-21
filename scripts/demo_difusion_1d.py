"""Ejecuta y grafica un caso configurable de difusion termica 1D.

Edita la seccion ``PARAMETROS EDITABLES`` y vuelve a ejecutar el archivo.
"""
import csv
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


RAIZ_PROYECTO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ_PROYECTO / "src"))

from mcplan.thermal.diffusion_1d import (  # noqa: E402
    calcular_difusividad_termica,
    simular_difusion_1d,
)


# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES: cambia estos valores para definir otro caso.
# Unidades: metros, segundos, °C, W/(m K), kg/m³ y J/(kg K).
# ---------------------------------------------------------------------------
LONGITUD_M = 1.0
NUMERO_NODOS = 41

TEMPERATURA_INICIAL_CENTRO_C = 0.0
TEMPERATURA_INICIAL_BASE_C = 0.0

# En cada extremo selecciona "dirichlet", "neumann" o "robin".
# Dirichlet: el valor asociado es temperatura fija [°C].
# Neumann: el valor asociado es gradiente dT/dx [K/m]; 0 significa aislado.
# Robin usa la temperatura ambiente definida más abajo y requiere h.

CONDICION_IZQUIERDA = "dirichlet"
VALOR_IZQUIERDO = 0.0
CONDICION_DERECHA = "dirichlet"
VALOR_DERECHO = 0.0
COEF_CONVECCION_IZQUIERDO_W_M2_K = 10.0
COEF_CONVECCION_DERECHO_W_M2_K = 10.0

# Ambiente sinusoidal configurable independientemente en cada lado.
# El lado izquierdo usa 20 ± 10 °C con período de 12 horas.
TEMPERATURA_MEDIA_AIRE_IZQUIERDO_C = 20.0
AMPLITUD_AIRE_IZQUIERDO_C = 10.0
PERIODO_AIRE_IZQUIERDO_S = 12.0 * 3600.0

# El lado derecho queda constante a 20 °C mientras la amplitud sea cero.
TEMPERATURA_MEDIA_AIRE_DERECHO_C = 20.0
AMPLITUD_AIRE_DERECHO_C = 0.0
PERIODO_AIRE_DERECHO_S = 12.0 * 3600.0


def temperatura_ambiente_izquierda_c(tiempo_s: float) -> float:
    """Temperatura del aire izquierda; varía de forma sinusoidal en el tiempo."""
    return TEMPERATURA_MEDIA_AIRE_IZQUIERDO_C + AMPLITUD_AIRE_IZQUIERDO_C * math.sin(
        2.0 * math.pi * tiempo_s / PERIODO_AIRE_IZQUIERDO_S
    )


def temperatura_ambiente_derecha_c(tiempo_s: float) -> float:
    """Temperatura del aire derecha; amplitud cero la mantiene constante."""
    return TEMPERATURA_MEDIA_AIRE_DERECHO_C + AMPLITUD_AIRE_DERECHO_C * math.sin(
        2.0 * math.pi * tiempo_s / PERIODO_AIRE_DERECHO_S
    )

# El solver calcula alfa = k / (rho * cp). Estos son valores ilustrativos.
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0

# Generación interna uniforme Q0, por ejemplo calor de hidratación [W/m³].
# Déjala en cero para simular difusión sin generación de calor.
GENERACION_CALOR_W_M3 = 500.0

TIEMPO_FINAL_S = 20000.0
NUMERO_PERFILES_GRAFICADOS = 8
TIEMPOS_GRAFICADOS_S = tuple(
    np.linspace(0.0, TIEMPO_FINAL_S, NUMERO_PERFILES_GRAFICADOS)
)
NUMERO_FOURIER_MAXIMO = 0.45  # Para FTCS debe ser <= 0.5.
# ---------------------------------------------------------------------------


def crear_perfil_inicial(posicion_m: np.ndarray) -> np.ndarray:
    """Crea un perfil seno sobre una temperatura base configurable."""
    fraccion = posicion_m / LONGITUD_M
    amplitud = TEMPERATURA_INICIAL_CENTRO_C - TEMPERATURA_INICIAL_BASE_C
    return TEMPERATURA_INICIAL_BASE_C + amplitud * np.sin(np.pi * fraccion)


def solucion_analitica_seno(
    posicion_m: np.ndarray, tiempo_s: float, difusividad_m2_s: float
) -> np.ndarray:
    """Solucion exacta del modo seno para extremos iguales y fijos."""
    longitud = LONGITUD_M
    amplitud = TEMPERATURA_INICIAL_CENTRO_C - TEMPERATURA_INICIAL_BASE_C
    return TEMPERATURA_INICIAL_BASE_C + amplitud * np.sin(
        np.pi * posicion_m / longitud
    ) * np.exp(-difusividad_m2_s * (np.pi / longitud) ** 2 * tiempo_s)


if NUMERO_NODOS < 3:
    raise ValueError("NUMERO_NODOS debe ser al menos 3")
if LONGITUD_M <= 0 or TIEMPO_FINAL_S <= 0:
    raise ValueError("LONGITUD_M y TIEMPO_FINAL_S deben ser positivos")
if not TIEMPOS_GRAFICADOS_S or any(
    t < 0 or t > TIEMPO_FINAL_S for t in TIEMPOS_GRAFICADOS_S
):
    raise ValueError("Cada tiempo graficado debe estar entre 0 y TIEMPO_FINAL_S")
for condicion in (CONDICION_IZQUIERDA, CONDICION_DERECHA):
    if condicion not in {"dirichlet", "neumann", "robin"}:
        raise ValueError("Las condiciones deben ser 'dirichlet', 'neumann' o 'robin'")

comparar_analitica = (
    GENERACION_CALOR_W_M3 == 0.0
    and CONDICION_IZQUIERDA == "dirichlet"
    and CONDICION_DERECHA == "dirichlet"
    and np.isclose(VALOR_IZQUIERDO, VALOR_DERECHO)
    and np.isclose(TEMPERATURA_INICIAL_BASE_C, VALOR_IZQUIERDO)
)

DIFUSIVIDAD_M2_S = calcular_difusividad_termica(
    CONDUCTIVIDAD_W_M_K,
    DENSIDAD_KG_M3,
    CALOR_ESPECIFICO_J_KG_K,
)
POSICION_M = np.linspace(0.0, LONGITUD_M, NUMERO_NODOS)
TEMPERATURA_INICIAL_C = crear_perfil_inicial(POSICION_M)
PASO_ESPACIAL_M = POSICION_M[1] - POSICION_M[0]
factores_estabilidad = [NUMERO_FOURIER_MAXIMO]
for condicion, coeficiente_h in (
    (CONDICION_IZQUIERDA, COEF_CONVECCION_IZQUIERDO_W_M2_K),
    (CONDICION_DERECHA, COEF_CONVECCION_DERECHO_W_M2_K),
):
    if condicion == "robin":
        biot_malla = coeficiente_h * PASO_ESPACIAL_M / CONDUCTIVIDAD_W_M_K
        factores_estabilidad.append(0.9 * 0.5 / (1.0 + biot_malla))
PASO_TIEMPO_MAX_S = (
    min(factores_estabilidad) * PASO_ESPACIAL_M**2 / DIFUSIVIDAD_M2_S
)

CARPETA_RESULTADOS = RAIZ_PROYECTO / "results" / "solver_1d"
CARPETA_RESULTADOS.mkdir(parents=True, exist_ok=True)
RUTA_GRAFICO = CARPETA_RESULTADOS / "perfiles_temperatura_1d.png"
RUTA_DATOS = CARPETA_RESULTADOS / "temperaturas_1d.csv"

colores = plt.cm.viridis(np.linspace(0.08, 0.88, len(TIEMPOS_GRAFICADOS_S)))
filas_csv = []
figura, ejes = plt.subplots(figsize=(9, 5.5))

for tiempo_s, color in zip(TIEMPOS_GRAFICADOS_S, colores):
    if tiempo_s == 0.0:
        # Mostrar la condición inicial pura, antes de imponer los bordes.
        # Así T(x, 0) puede diferir del borde Dirichlet para t > 0.
        temperatura_c = TEMPERATURA_INICIAL_C.copy()
    else:
        resultado = simular_difusion_1d(
            temperatura_inicial_c=TEMPERATURA_INICIAL_C,
            longitud_m=LONGITUD_M,
            difusividad_m2_s=DIFUSIVIDAD_M2_S,
            condicion_izquierda=CONDICION_IZQUIERDA,
            valor_izquierdo=(
                VALOR_IZQUIERDO if CONDICION_IZQUIERDA != "robin" else None
            ),
            condicion_derecha=CONDICION_DERECHA,
            valor_derecho=(
                VALOR_DERECHO if CONDICION_DERECHA != "robin" else None
            ),
            conductividad_w_m_k=CONDUCTIVIDAD_W_M_K,
            coeficiente_conveccion_izquierdo_w_m2_k=(
                COEF_CONVECCION_IZQUIERDO_W_M2_K
            ),
            coeficiente_conveccion_derecho_w_m2_k=(
                COEF_CONVECCION_DERECHO_W_M2_K
            ),
            temperatura_ambiente_izquierda_c=temperatura_ambiente_izquierda_c,
            temperatura_ambiente_derecha_c=temperatura_ambiente_derecha_c,
            generacion_calor_w_m3=GENERACION_CALOR_W_M3,
            capacidad_volumetrica_j_m3_k=(
                DENSIDAD_KG_M3 * CALOR_ESPECIFICO_J_KG_K
            ),
            tiempo_final_s=tiempo_s,
            paso_tiempo_max_s=PASO_TIEMPO_MAX_S,
        )
        temperatura_c = resultado.temperatura_c

    ejes.plot(
        POSICION_M,
        temperatura_c,
        color=color,
        linewidth=2.0,
        marker="o",
        markersize=3.0,
        markevery=2,
        label=f"Numérica, t = {tiempo_s / 60.0:.0f} min",
    )

    if comparar_analitica:
        exacta = solucion_analitica_seno(
            POSICION_M, tiempo_s, DIFUSIVIDAD_M2_S
        )
        if tiempo_s == TIEMPOS_GRAFICADOS_S[-1]:
            ejes.plot(
                POSICION_M,
                exacta,
                color="#24384B",
                linestyle="--",
                linewidth=1.5,
                label="Analítica al final",
            )
    else:
        exacta = np.full_like(POSICION_M, np.nan)

    filas_csv.extend(
        (tiempo_s, posicion, temperatura_num, temperatura_exacta)
        for posicion, temperatura_num, temperatura_exacta
        in zip(POSICION_M, temperatura_c, exacta)
    )

error_maximo_c = np.nanmax(
    np.abs(np.asarray([fila[2] - fila[3] for fila in filas_csv]))
) if comparar_analitica else None

if "robin" in (CONDICION_IZQUIERDA, CONDICION_DERECHA):
    titulo = "Difusión térmica 1D · Robin con ambiente variable"
else:
    titulo = "Difusión térmica 1D · caso configurable"
ejes.set_title(titulo)
ejes.set_xlabel("Posición en la barra, x [m]")
ejes.set_ylabel("Temperatura [°C]")
ejes.grid(True, alpha=0.28)
ejes.legend(ncol=2, fontsize=8.5, frameon=False)
figura.tight_layout()
figura.savefig(RUTA_GRAFICO, dpi=180, bbox_inches="tight")

with RUTA_DATOS.open("w", newline="", encoding="utf-8") as archivo:
    escritor = csv.writer(archivo)
    escritor.writerow(
        ("tiempo_s", "posicion_m", "temperatura_numerica_c", "temperatura_analitica_c")
    )
    escritor.writerows(filas_csv)

print(f"Difusividad calculada: {DIFUSIVIDAD_M2_S:.3e} m^2/s")
print(f"Generación interna Q0: {GENERACION_CALOR_W_M3:.3g} W/m^3")
print(f"Temperatura inicial en el centro: {TEMPERATURA_INICIAL_CENTRO_C:.2f} °C")
print(
    f"Perfil inicial: base {TEMPERATURA_INICIAL_BASE_C:.2f} °C, "
    f"centro {TEMPERATURA_INICIAL_CENTRO_C:.2f} °C"
)
if CONDICION_IZQUIERDA == "robin":
    print("Borde izquierdo: Robin")
else:
    print(f"Borde izquierdo: {CONDICION_IZQUIERDA} = {VALOR_IZQUIERDO}")
if CONDICION_DERECHA == "robin":
    print("Borde derecho: Robin")
else:
    print(f"Borde derecho: {CONDICION_DERECHA} = {VALOR_DERECHO}")
if CONDICION_IZQUIERDA == "robin":
    print(
        "Aire izquierdo: "
        f"{temperatura_ambiente_izquierda_c(0.0):.2f} °C al inicio, "
        f"{temperatura_ambiente_izquierda_c(TIEMPO_FINAL_S):.2f} °C al final"
    )
    print(f"h izquierdo: {COEF_CONVECCION_IZQUIERDO_W_M2_K:.3g} W/(m^2 K)")
if CONDICION_DERECHA == "robin":
    print(
        "Aire derecho: "
        f"{temperatura_ambiente_derecha_c(0.0):.2f} °C al inicio, "
        f"{temperatura_ambiente_derecha_c(TIEMPO_FINAL_S):.2f} °C al final"
    )
    print(f"h derecho: {COEF_CONVECCION_DERECHO_W_M2_K:.3g} W/(m^2 K)")
print(f"Nodos: {NUMERO_NODOS}; paso espacial: {PASO_ESPACIAL_M:.4f} m")
if error_maximo_c is not None:
    print(f"Error máximo contra la solución analítica: {error_maximo_c:.6f} °C")
else:
    print("Comparación analítica omitida: los datos ya no cumplen el caso seno de referencia.")
print(f"Gráfico guardado en: {RUTA_GRAFICO}")
print(f"Datos guardados en: {RUTA_DATOS}")

plt.show()
plt.close(figura)
