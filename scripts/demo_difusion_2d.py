r"""Demo de difusion termica 2D (solver explicito por balance de celda).

Edita ``PARAMETROS EDITABLES`` y ejecuta:  python .\scripts\demo_difusion_2d.py
Parte 1: caso base 2D (por defecto, flujo esencialmente unidireccional en x); sus
caras se eligen en BORDES DEL CASO BASE. Compara escenarios que cambian UNA cosa
a la vez: dx, dt, una condicion de borde y la fuente. Es un bloque rectangular de prueba (NO la fundacion
voxelizada). Unidades: m, s, °C, W/(m K), kg/m³, J/(kg K), W/m³, W/m².
"""
import csv
from pathlib import Path
import sys

# True: abre ventanas de Matplotlib al final (una sola llamada a plt.show()).
# False: sin ventanas (backend Agg), para ejecuciones automaticas; igual guarda los PNG.
MOSTRAR_GRAFICOS = True

import matplotlib
if not MOSTRAR_GRAFICOS:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RAIZ_PROYECTO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ_PROYECTO / "src"))

from mcplan.thermal.diffusion_nd import (  # noqa: E402
    DifusionCartesiana, aislado, dirichlet, interpolar_en_punto, neumann, robin)

# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES
# ---------------------------------------------------------------------------
LADO_X_M, LADO_Y_M = 0.50, 0.30        # tamano del bloque [m]
CELDAS_X, CELDAS_Y = 25, 15            # celdas base (dx = LADO/CELDAS)
PROFUNDIDAD_M = 1.0                    # solo 2D: espesor fuera del plano [m]
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0
TEMPERATURA_INICIAL_C = 20.0
TIEMPO_FINAL_S = 6 * 3600.0
FRACCION_PASO_ESTABLE = 0.8            # dt = fraccion * limite de estabilidad
FUENTE_W_M3 = 0.0                      # escenario base [W/m³]
FUENTE_ESCENARIO_W_M3 = 200.0          # escenario "con fuente"
CARPETA_SALIDA = RAIZ_PROYECTO / "results" / "solver_nd"

# ---- BORDES DEL CASO BASE (Parte 1: escenarios dx, dt, borde y fuente) ----
# Una variable por cara. MIN/MAX = cara de coordenada minima/maxima de ese eje.
# Cada una acepta directamente:
#   dirichlet(20.0)   temperatura fija de 20 °C en la cara
#   neumann(100.0)    flujo ENTRANTE de 100 W/m²  (positivo calienta)
#   neumann(-100.0)   flujo SALIENTE de 100 W/m²  (negativo enfria)
#   aislado()         flujo cero (adiabatica)
#   robin(10.0, 5.0)  conveccion: h = 10 W/(m² K), ambiente a 5 °C
BORDE_X_MIN = dirichlet(60.0)          # x = 0
BORDE_X_MAX = robin(10.0, 5.0)         # x = LADO_X_M  (h [W/(m² K)], T ambiente [°C])
BORDE_Y_MIN = aislado()                # y = 0
BORDE_Y_MAX = aislado()                # y = LADO_Y_M
# El escenario "x+ aislado" REEMPLAZA a proposito la cara BORDE_X_MAX del caso
# base por aislado(); las demas caras se mantienen como las definiste arriba.

# ---- BORDES DEL CASO MULTIDIMENSIONAL (Parte 2), independientes del caso base ----
# Mismas opciones y unidades que arriba.
MD_BORDE_X_MIN = dirichlet(50.0)       # x = 0
MD_BORDE_X_MAX = robin(10.0, 5.0)      # x = LADO_X_M
MD_BORDE_Y_MIN = dirichlet(10.0)       # y = 0
MD_BORDE_Y_MAX = neumann(100.0)        # y = LADO_Y_M: flujo ENTRANTE de 100 W/m²
# Temperatura inicial localizada (mancha gaussiana) del caso multidimensional.
MANCHA_M = (0.15, 0.10)              # centro de la mancha [m] (x, y)
SIGMA_MANCHA_M = 0.05                  # ancho de la mancha inicial [m]
T_MANCHA_C = 80.0                      # temperatura en el centro de la mancha [°C]
UMBRAL_VARIACION_C = 1.0               # variacion minima para considerar que T varia en un eje [K]
# El paso de tiempo estable se calcula solo a partir de las caras elegidas.


def fronteras_base():
    """Caras del caso base, tomadas de las variables BORDE_* de PARAMETROS EDITABLES."""
    return {"x-": BORDE_X_MIN, "x+": BORDE_X_MAX, "y-": BORDE_Y_MIN, "y+": BORDE_Y_MAX}


def correr(nombre, *, factor_celdas=1, factor_dt=1.0, fronteras=None, fuente=FUENTE_W_M3):
    nx, ny = CELDAS_X * factor_celdas, CELDAS_Y * factor_celdas
    sistema = DifusionCartesiana(
        forma=(nx, ny), dx_m=(LADO_X_M / nx, LADO_Y_M / ny),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K, densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K, fuente_w_m3=fuente,
        fronteras=fronteras or fronteras_base(), profundidad_m=PROFUNDIDAD_M)
    paso = FRACCION_PASO_ESTABLE * sistema.paso_maximo_estable_s() * factor_dt
    r = sistema.simular(temperatura_inicial_c=TEMPERATURA_INICIAL_C,
                        tiempo_final_s=TIEMPO_FINAL_S, paso_tiempo_max_s=paso)
    T = r.temperatura_c
    return {"escenario": nombre, "dx_m": LADO_X_M / nx, "dt_s": r.paso_tiempo_s,
            "dt_limite_s": r.paso_maximo_estable_s, "pasos": r.numero_pasos,
            "T_centro_C": interpolar_en_punto(T, (LADO_X_M / nx, LADO_Y_M / ny), (LADO_X_M / 2, LADO_Y_M / 2)), "T_min_C": float(T.min()),
            "T_max_C": float(T.max()), "error_balance_J": r.error_balance_j,
            "campo": T}


def main():
    escenarios = [
        correr("base"),
        correr("dx/2", factor_celdas=2),
        correr("dt/2", factor_dt=0.5),
        correr("x+ aislado", fronteras={**fronteras_base(), "x+": aislado()}),
        correr("con fuente", fuente=FUENTE_ESCENARIO_W_M3),
    ]
    print(f"2D {CELDAS_X}x{CELDAS_Y} celdas base, t = {TIEMPO_FINAL_S/3600:g} h, "
          f"T0 = {TEMPERATURA_INICIAL_C} °C (bloque de prueba, no la fundacion)")
    print(f"{'escenario':<12}{'dx[m]':>8}{'dt[s]':>9}{'dt_lim[s]':>11}{'pasos':>7}"
          f"{'T(centro fis.)':>15}{'Tmin':>8}{'Tmax':>8}{'err.balance[J]':>16}")
    for e in escenarios:
        print(f"{e['escenario']:<12}{e['dx_m']:>8.4f}{e['dt_s']:>9.2f}{e['dt_limite_s']:>11.2f}"
              f"{e['pasos']:>7d}{e['T_centro_C']:>15.3f}{e['T_min_C']:>8.2f}"
              f"{e['T_max_C']:>8.2f}{e['error_balance_J']:>16.3e}")

    print("T(centro fis.) = interpolacion multilineal en el centro geometrico del bloque (mismo punto en todas las mallas).")
    print("Nota: 'dx/2' tambien reduce dt (~dx^2) por estabilidad: es refinamiento CONJUNTO, no espacial puro.")
    print("Nota (caso base): " + describir_direccionalidad(escenarios[0]["campo"]))
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)
    with open(CARPETA_SALIDA / "resumen_2d.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        claves = [k for k in escenarios[0] if k != "campo"]
        w.writerow(claves)
        for e in escenarios:
            w.writerow([e[k] for k in claves])
    fig, ejes = plt.subplots(1, len(escenarios), figsize=(3.2 * len(escenarios), 2.8),
                             sharey=True, constrained_layout=True)
    vmin = min(e["T_min_C"] for e in escenarios)
    vmax = max(e["T_max_C"] for e in escenarios)
    for eje, e in zip(ejes, escenarios):
        im = eje.imshow(e["campo"].T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
                        vmin=vmin, vmax=vmax, cmap="inferno")
        eje.set_title(e["escenario"]); eje.set_xlabel("x [m]")
    ejes[0].set_ylabel("y [m]")
    fig.colorbar(im, ax=ejes, label="T [°C]")
    fig.savefig(CARPETA_SALIDA / "campos_2d.png", dpi=110)
    print(f"Imagen guardada: {CARPETA_SALIDA / 'campos_2d.png'}")
    reportar_multidimensional()

# ---------------------------------------------------------------------------
# PARTE 2: caso claramente MULTIDIMENSIONAL (gradientes en TODOS los ejes)
# Temperatura inicial localizada (mancha caliente) + caras con condiciones
# distintas. Se comprueba numericamente que T varia en cada direccion.
# ---------------------------------------------------------------------------


def variacion_por_eje(T):
    """Mayor diferencia max-min de T a lo largo de cada eje (otras posiciones fijas)."""
    return [float(np.max(np.ptp(T, axis=a))) for a in range(T.ndim)]


def describir_direccionalidad(T):
    """Informa (sin abortar) en cuantos ejes T varia mas de UMBRAL_VARIACION_C."""
    ejes_txt = "xyz"
    sin = [ejes_txt[a] for a, v in enumerate(variacion_por_eje(T)) if v <= UMBRAL_VARIACION_C]
    if not sin:
        return "T varia en todos los ejes: caso multidimensional."
    return ("sin variacion apreciable (> %g K) en: %s. Con esta configuracion el caso es "
            "uniforme o unidireccional; se guardan y muestran los graficos igualmente."
            % (UMBRAL_VARIACION_C, ", ".join(sin)))


def caso_multidimensional():
    nx, ny = CELDAS_X, CELDAS_Y
    sistema = DifusionCartesiana(
        forma=(nx, ny), dx_m=(LADO_X_M / nx, LADO_Y_M / ny),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K, densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K, profundidad_m=PROFUNDIDAD_M,
        fronteras={"x-": MD_BORDE_X_MIN, "x+": MD_BORDE_X_MAX,
                   "y-": MD_BORDE_Y_MIN, "y+": MD_BORDE_Y_MAX})
    x = (np.arange(nx) + 0.5) * LADO_X_M / nx
    y = (np.arange(ny) + 0.5) * LADO_Y_M / ny
    X, Y = np.meshgrid(x, y, indexing="ij")
    T0 = TEMPERATURA_INICIAL_C + (T_MANCHA_C - TEMPERATURA_INICIAL_C) * np.exp(
        -((X - MANCHA_M[0]) ** 2 + (Y - MANCHA_M[1]) ** 2) / (2 * SIGMA_MANCHA_M ** 2))
    r = sistema.simular(temperatura_inicial_c=T0, tiempo_final_s=TIEMPO_FINAL_S,
                        paso_tiempo_max_s=FRACCION_PASO_ESTABLE * sistema.paso_maximo_estable_s())
    return sistema, r, T0


def reportar_multidimensional():
    sistema, r, T0 = caso_multidimensional()
    T = r.temperatura_c
    ejes_txt = "xyz"
    print("\nCaso multidimensional: mancha caliente + caras distintas "
          f"(dt = {r.paso_tiempo_s:.2f} s, limite {r.paso_maximo_estable_s:.2f} s, "
          f"err. balance {r.error_balance_j:.2e} J)")
    for a, v in enumerate(variacion_por_eje(T)):
        estado = "variacion" if v > UMBRAL_VARIACION_C else "sin variacion apreciable"
        print(f"  variacion de T en {ejes_txt[a]}: {v:8.2f} K  [{estado}]")
    print("  Nota: " + describir_direccionalidad(T))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3), constrained_layout=True)
    for eje, campo, titulo in ((a1, T0, "inicial"), (a2, T, f"t = {TIEMPO_FINAL_S/3600:g} h")):
        im = eje.imshow(campo.T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
                        vmin=min(T0.min(), T.min()), vmax=max(T0.max(), T.max()), cmap="inferno")
        eje.set_title(titulo); eje.set_xlabel("x [m]"); eje.set_ylabel("y [m]")
    fig.colorbar(im, ax=[a1, a2], label="T [°C]")
    fig.savefig(CARPETA_SALIDA / "multidimensional_2d.png", dpi=110)
    print(f"Imagen guardada: {CARPETA_SALIDA / 'multidimensional_2d.png'}")


if __name__ == "__main__":
    main()
    if MOSTRAR_GRAFICOS:
        plt.show()   # una sola llamada: muestra todas las figuras ya guardadas
