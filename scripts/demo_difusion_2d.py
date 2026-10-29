r"""Demo 2D: bloque 4 x 2 m con generación interna Qdot(t).

Ejecuta: python .\scripts\demo_difusion_2d.py
Temperatura inicial y las cuatro caras: 20 °C. Qdot(t) sube hasta Qmax en T=1 día,
baja a cero en 2T y se mantiene en cero después. Guarda mapas, historia T(t),
datos CSV y una animación GIF.
Es un bloque rectangular de prueba, no la fundación voxelizada. Unidades SI.
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
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

RAIZ_PROYECTO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ_PROYECTO / "src"))

from mcplan.thermal.diffusion_nd import (  # noqa: E402
    DifusionCartesiana, aislado, dirichlet, interpolar_en_punto, neumann, robin)

# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES: caso solicitado (bloque 4 x 2 m)
# ---------------------------------------------------------------------------
LADO_X_M, LADO_Y_M = 4.0, 2.0          # tamano del bloque [m]
CELDAS_X, CELDAS_Y = 80, 40            # dx = dy = 0.05 m
PROFUNDIDAD_M = 1.0                    # solo 2D: espesor fuera del plano [m]
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0
TEMPERATURA_INICIAL_C = 20.0
TIEMPO_FINAL_S = 50 * 3600.0       # observar pulso y enfriamiento posterior
FRACCION_PASO_ESTABLE = 0.8            # dt = fraccion * limite de estabilidad
QMAX_W_M3 = 1300.0                      # valor ilustrativo y editable [W/m³]
TIEMPO_PEAK_Q_S = 24.0 * 3600.0        # T de la función: el peak ocurre al día 1
DURACION_FUENTE_S = 2.0 * TIEMPO_PEAK_Q_S  # Q(t) vuelve a cero al día 2
NUMERO_INTERVALOS_ANIMACION = 25      # un cuadro por hora durante los 10 días
VELOCIDAD_GIF_FPS = 8
CARPETA_SALIDA = RAIZ_PROYECTO / "results" / "solver_nd"

# ---- BORDES DEL CASO SOLICITADO ----
# Una variable por cara. MIN/MAX = cara de coordenada minima/maxima de ese eje.
# Cada una acepta directamente:
#   dirichlet(20.0)   temperatura fija de 20 °C en la cara
#   neumann(100.0)    flujo ENTRANTE de 100 W/m²  (positivo calienta)
#   neumann(-100.0)   flujo SALIENTE de 100 W/m²  (negativo enfria)
#   aislado()         flujo cero (adiabatica)
#   robin(10.0, 5.0)  conveccion: h = 10 W/(m² K), ambiente a 5 °C
BORDE_X_MIN = dirichlet(20.0)          # x = 0
BORDE_X_MAX = dirichlet(20.0)          # x = LADO_X_M
BORDE_Y_MIN = dirichlet(20.0)          # y = 0
BORDE_Y_MAX = dirichlet(20.0)          # y = LADO_Y_M
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


def fuente_calor_w_m3(tiempo_s: float) -> float:
    """Q(t)=Qmax/2*(1-cos(pi*t/T)) hasta 2T; después queda en cero."""
    if tiempo_s < 0.0 or tiempo_s > DURACION_FUENTE_S:
        return 0.0
    return (QMAX_W_M3 / 2.0) * (
        1.0 - np.cos(np.pi * tiempo_s / TIEMPO_PEAK_Q_S)
    )


def crear_sistema_base():
    return DifusionCartesiana(
        forma=(CELDAS_X, CELDAS_Y),
        dx_m=(LADO_X_M / CELDAS_X, LADO_Y_M / CELDAS_Y),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K,
        densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K,
        fuente_w_m3=fuente_calor_w_m3,
        fronteras=fronteras_base(),
        profundidad_m=PROFUNDIDAD_M,
    )


def simular_historial():
    """Devuelve el campo térmico y temperaturas representativas en cada cuadro."""
    sistema = crear_sistema_base()
    tiempos_s = np.linspace(0.0, TIEMPO_FINAL_S, NUMERO_INTERVALOS_ANIMACION + 1)
    temperatura = np.full((CELDAS_X, CELDAS_Y), TEMPERATURA_INICIAL_C, dtype=float)
    campos = [temperatura.copy()]
    paso_max_s = FRACCION_PASO_ESTABLE * sistema.paso_maximo_estable_s()
    errores_balance = []

    for inicio_s, fin_s in zip(tiempos_s[:-1], tiempos_s[1:]):
        resultado = sistema.simular(
            temperatura_inicial_c=temperatura,
            tiempo_final_s=float(fin_s - inicio_s),
            paso_tiempo_max_s=paso_max_s,
            tiempo_inicial_s=float(inicio_s),
        )
        temperatura = resultado.temperatura_c
        campos.append(temperatura.copy())
        errores_balance.append(resultado.error_balance_j)

    x_centro_m = (LADO_X_M / 2.0, LADO_Y_M / 2.0)
    dx_m = (LADO_X_M / CELDAS_X, LADO_Y_M / CELDAS_Y)
    historia = {
        "tiempos_s": tiempos_s,
        "campos": campos,
        "centro_c": np.array([interpolar_en_punto(T, dx_m, x_centro_m) for T in campos]),
        "promedio_c": np.array([float(T.mean()) for T in campos]),
        "minimo_c": np.array([float(T.min()) for T in campos]),
        "maximo_c": np.array([float(T.max()) for T in campos]),
        "fuente_w_m3": np.array([fuente_calor_w_m3(float(t)) for t in tiempos_s]),
        "dt_limite_s": sistema.paso_maximo_estable_s(),
        "dt_max_s": paso_max_s,
        "error_balance_j": float(sum(errores_balance)),
    }
    return historia


def guardar_mapa(campo, tiempo_s, ruta, vmin, vmax):
    fig, eje = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    imagen = eje.imshow(
        campo.T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
        vmin=vmin, vmax=vmax, cmap="inferno", aspect="auto",
    )
    eje.set_xlabel("x [m]")
    eje.set_ylabel("y [m]")
    eje.set_title(f"Temperatura del bloque · t = {tiempo_s / 3600:g} h")
    fig.colorbar(imagen, ax=eje, label="Temperatura [°C]")
    fig.savefig(ruta, dpi=160)
    return fig


def guardar_historia_estatica(historia, ruta):
    horas = historia["tiempos_s"] / 3600.0
    indice_max = int(np.argmax(historia["centro_c"]))
    fig, eje = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    eje.plot(horas, historia["centro_c"], marker="o", ms=3,
             label="Temperatura en el centro")
    eje.plot(horas, historia["promedio_c"], marker="s", ms=3,
             label="Temperatura promedio del bloque")
    eje.axvline(TIEMPO_PEAK_Q_S / 3600.0, color="tab:red", linestyle=":",
                alpha=0.75, label="Qmax · día 1")
    eje.axvline(DURACION_FUENTE_S / 3600.0, color="black", linestyle="--",
                alpha=0.65, label="Q(t) se apaga · día 2")
    eje.scatter(horas[indice_max], historia["centro_c"][indice_max], marker="*",
                s=110, color="tab:blue", zorder=5, label="Máximo del centro")
    eje.set_xlabel("Tiempo [h]")
    eje.set_ylabel("Temperatura [°C]")
    eje.set_title("Evolución de la temperatura")
    eje.grid(True, alpha=0.35)
    eje.legend()
    fig.savefig(ruta, dpi=160)
    return fig


def guardar_animacion(historia, ruta, vmin, vmax):
    tiempos_h = historia["tiempos_s"] / 3600.0
    campos = historia["campos"]
    fig, (eje_mapa, eje_historia) = plt.subplots(
        1, 2, figsize=(11, 4.6), constrained_layout=True,
    )
    imagen = eje_mapa.imshow(
        campos[0].T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
        vmin=vmin, vmax=vmax, cmap="inferno", aspect="auto",
        animated=True,
    )
    eje_mapa.set_xlabel("x [m]")
    eje_mapa.set_ylabel("y [m]")
    titulo_mapa = eje_mapa.set_title("Temperatura · t = 0 h")
    fig.colorbar(imagen, ax=eje_mapa, label="Temperatura [°C]")

    eje_historia.plot(tiempos_h, historia["centro_c"], color="tab:blue",
                      alpha=0.35, label="Centro")
    eje_historia.plot(tiempos_h, historia["promedio_c"], color="tab:orange",
                      alpha=0.35, label="Promedio")
    eje_historia.axvline(TIEMPO_PEAK_Q_S / 3600.0, color="tab:red",
                         linestyle=":", alpha=0.75, label="Qmax · día 1")
    eje_historia.axvline(DURACION_FUENTE_S / 3600.0, color="black",
                         linestyle="--", alpha=0.65, label="Q(t) se apaga · día 2")
    linea_centro, = eje_historia.plot([], [], color="tab:blue", lw=2,
                                      label="Centro (avance)")
    linea_promedio, = eje_historia.plot([], [], color="tab:orange", lw=2,
                                        label="Promedio (avance)")
    punto_centro, = eje_historia.plot([], [], marker="o", color="tab:blue")
    punto_promedio, = eje_historia.plot([], [], marker="o", color="tab:orange")
    eje_historia.set_xlim(0.0, TIEMPO_FINAL_S / 3600.0)
    eje_historia.set_ylim(vmin - 0.5, vmax + 0.5)
    eje_historia.set_xlabel("Tiempo [h]")
    eje_historia.set_ylabel("Temperatura [°C]")
    eje_historia.set_title("Temperatura vs. tiempo")
    eje_historia.grid(True, alpha=0.35)
    eje_historia.legend(fontsize=8)

    def actualizar(indice):
        imagen.set_data(campos[indice].T)
        titulo_mapa.set_text(f"Temperatura · t = {tiempos_h[indice]:.2f} h")
        linea_centro.set_data(tiempos_h[:indice + 1], historia["centro_c"][:indice + 1])
        linea_promedio.set_data(tiempos_h[:indice + 1], historia["promedio_c"][:indice + 1])
        punto_centro.set_data([tiempos_h[indice]], [historia["centro_c"][indice]])
        punto_promedio.set_data([tiempos_h[indice]], [historia["promedio_c"][indice]])
        return (imagen, titulo_mapa, linea_centro, linea_promedio,
                punto_centro, punto_promedio)

    pelicula = FuncAnimation(fig, actualizar, frames=len(campos), interval=350,
                              blit=False, repeat=True)
    pelicula.save(ruta, writer=PillowWriter(fps=VELOCIDAD_GIF_FPS))
    plt.close(fig)


def main():
    CARPETA_SALIDA.mkdir(parents=True, exist_ok=True)
    historia = simular_historial()
    campos = historia["campos"]
    tiempos_h = historia["tiempos_s"] / 3600.0
    vmin = min(float(T.min()) for T in campos)
    vmax = max(float(T.max()) for T in campos)

    guardar_mapa(
        campos[0], 0.0, CARPETA_SALIDA / "temperatura_2d_t0.png", vmin, vmax,
    )
    indice_6h = int(np.argmin(np.abs(historia["tiempos_s"] - 6.0 * 3600.0)))
    guardar_mapa(
        campos[indice_6h], historia["tiempos_s"][indice_6h],
        CARPETA_SALIDA / "temperatura_2d_t6h.png", vmin, vmax,
    )
    dias_finales = TIEMPO_FINAL_S / (24.0 * 3600.0)
    ruta_final = CARPETA_SALIDA / f"temperatura_2d_t{dias_finales:g}d.png"
    final = guardar_mapa(campos[-1], TIEMPO_FINAL_S, ruta_final, vmin, vmax)
    # Mantiene el nombre histórico del mapa final para compatibilidad.
    final.savefig(CARPETA_SALIDA / "campos_2d.png", dpi=160)
    guardar_historia_estatica(
        historia, CARPETA_SALIDA / "historia_temperatura_2d.png",
    )
    guardar_animacion(
        historia, CARPETA_SALIDA / "evolucion_temperatura_2d.gif", vmin, vmax,
    )

    with open(CARPETA_SALIDA / "historia_temperatura_2d.csv", "w", newline="",
              encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(["tiempo_h", "Qdot_W_m3", "temperatura_centro_C",
                           "temperatura_promedio_C", "temperatura_minima_C",
                           "temperatura_maxima_C"])
        for i, hora in enumerate(tiempos_h):
            escritor.writerow([hora, historia["fuente_w_m3"][i], historia["centro_c"][i],
                               historia["promedio_c"][i], historia["minimo_c"][i],
                               historia["maximo_c"][i]])

    with open(CARPETA_SALIDA / "resumen_2d.csv", "w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        indice_max_centro = int(np.argmax(historia["centro_c"]))
        escritor.writerow(["LadoX_m", "LadoY_m", "CeldasX", "CeldasY", "Qmax_W_m3",
                           "TiempoPeakQ_h", "TiempoCorteQ_h", "TiempoFinal_h",
                           "TemperaturaCentroFinal_C",
                           "TemperaturaCentroMaxima_C", "TiempoCentroMaximo_h",
                           "TemperaturaPromedioFinal_C", "TemperaturaMinimaFinal_C",
                           "TemperaturaMaximaFinal_C", "PasoMaximo_s", "ErrorBalance_J"])
        escritor.writerow([LADO_X_M, LADO_Y_M, CELDAS_X, CELDAS_Y,
                           QMAX_W_M3, TIEMPO_PEAK_Q_S / 3600.0,
                           DURACION_FUENTE_S / 3600.0, TIEMPO_FINAL_S / 3600.0,
                           historia["centro_c"][-1], historia["centro_c"][indice_max_centro],
                           tiempos_h[indice_max_centro], historia["promedio_c"][-1],
                           historia["minimo_c"][-1], historia["maximo_c"][-1],
                           historia["dt_max_s"],
                           historia["error_balance_j"]])

    print(f"Bloque 2D: {LADO_X_M:g} x {LADO_Y_M:g} m; malla {CELDAS_X} x {CELDAS_Y}; "
          f"t final = {TIEMPO_FINAL_S / 86400:g} días; T inicial = {TEMPERATURA_INICIAL_C:g} °C")
    print("Bordes Dirichlet: 20 °C en x-, x+, y-, y+; fuente: semicoseno de la fórmula indicada")
    print(f"Qmax = {QMAX_W_M3:g} W/m³ a t = {TIEMPO_PEAK_Q_S / 3600:g} h; "
          f"Q(t) vuelve a cero a las {DURACION_FUENTE_S / 3600:g} h")
    print(f"Temperatura máxima del centro: {historia['centro_c'][indice_max_centro]:.3f} °C "
          f"a las {tiempos_h[indice_max_centro]:.2f} h")
    print(f"Temperatura en el centro al final: {historia['centro_c'][-1]:.3f} °C")
    print(f"Promedio al final: {historia['promedio_c'][-1]:.3f} °C; "
          f"mínimo/máximo: {historia['minimo_c'][-1]:.3f}/{historia['maximo_c'][-1]:.3f} °C")
    print(f"Paso estable máximo: {historia['dt_limite_s']:.2f} s; "
          f"error total de balance: {historia['error_balance_j']:.3e} J")
    print("Archivos guardados en:", CARPETA_SALIDA)

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
