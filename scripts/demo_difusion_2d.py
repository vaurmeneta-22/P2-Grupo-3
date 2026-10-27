r"""Demo de difusion termica 2D (solver explicito por balance de celda).

Edita ``PARAMETROS EDITABLES`` y ejecuta:  python .\scripts\demo_difusion_2d.py
Parte 1: flujo esencialmente UNIDIRECCIONAL (x) sobre una malla 2D. Compara escenarios que cambian UNA cosa a la vez: dx, dt, una condicion de
borde y la fuente. Es un bloque rectangular de prueba (NO la fundacion
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
H_CONVECCION_W_M2_K = 10.0             # cara x+ en el caso base (Robin)
T_AMBIENTE_C = 5.0
T_DIRICHLET_C = 60.0                   # cara x- en el caso base
CARPETA_SALIDA = RAIZ_PROYECTO / "results" / "solver_nd"


def fronteras_base():
    """x-: Dirichlet; x+: Robin; y-, y+: aisladas (Neumann q=0)."""
    return {"x-": dirichlet(T_DIRICHLET_C), "x+": robin(H_CONVECCION_W_M2_K, T_AMBIENTE_C),
            "y-": aislado(), "y+": aislado()}


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
    print("Nota: estos escenarios son casi unidireccionales (x): casi no hay variacion en y.")
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
MANCHA_M = (0.15, 0.10)              # centro de la mancha [m] (x, y)
SIGMA_MANCHA_M = 0.05                  # ancho de la mancha inicial [m]
T_MANCHA_C = 80.0                      # temperatura en el centro de la mancha
UMBRAL_VARIACION_C = 1.0               # variacion minima exigida por eje [K]


def variacion_por_eje(T):
    """Mayor diferencia max-min de T a lo largo de cada eje (otras posiciones fijas)."""
    return [float(np.max(np.ptp(T, axis=a))) for a in range(T.ndim)]


def caso_multidimensional():
    nx, ny = CELDAS_X, CELDAS_Y
    sistema = DifusionCartesiana(
        forma=(nx, ny), dx_m=(LADO_X_M / nx, LADO_Y_M / ny),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K, densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K, profundidad_m=PROFUNDIDAD_M,
        fronteras={"x-": dirichlet(50.0), "x+": robin(H_CONVECCION_W_M2_K, T_AMBIENTE_C),
                   "y-": dirichlet(10.0), "y+": neumann(100.0)})
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
        estado = "OK" if v > UMBRAL_VARIACION_C else "SIN VARIACION"
        print(f"  variacion de T en {ejes_txt[a]}: {v:8.2f} K  [{estado}]")
    assert all(v > UMBRAL_VARIACION_C for v in variacion_por_eje(T)), "falta variacion en algun eje"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3), constrained_layout=True)
    for eje, campo, titulo in ((a1, T0, "inicial"), (a2, T, f"t = {TIEMPO_FINAL_S/3600:g} h")):
        im = eje.imshow(campo.T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
                        vmin=10, vmax=T_MANCHA_C, cmap="inferno")
        eje.set_title(titulo); eje.set_xlabel("x [m]"); eje.set_ylabel("y [m]")
    fig.colorbar(im, ax=[a1, a2], label="T [°C]")
    fig.savefig(CARPETA_SALIDA / "multidimensional_2d.png", dpi=110)
    print(f"Imagen guardada: {CARPETA_SALIDA / 'multidimensional_2d.png'}")


if __name__ == "__main__":
    main()
    if MOSTRAR_GRAFICOS:
        plt.show()   # una sola llamada: muestra todas las figuras ya guardadas
