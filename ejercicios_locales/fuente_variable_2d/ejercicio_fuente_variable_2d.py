r"""Ejercicio de clases: fuente uniforme VARIABLE en el tiempo, rectangulo 2D.

Pulsa "Ejecutar archivo de Python" (VS Code) o:  python ejercicio_fuente_variable_2d.py
Funciona desde cualquier directorio. Es un experimento LOCAL, separado de E04/E05: usa su
propia copia del solver (solver_ejercicio.py, en esta carpeta).

Caso: rectangulo 4 m x 2 m, T inicial uniforme y cuatro bordes Dirichlet a 20 °C,
material homogeneo y fuente uniforme en el espacio
    Q(t) = Qmax/2 [1 - cos(2 pi t / tau)]   (0 <= t <= tau),    Q(t) = 0 fuera.
Objetivo: el mayor Qmax >= 0 tal que el centro NO supere T_LIMITE_C en todo el horizonte.
Como el problema es lineal en Qmax (T0 = temperatura de borde = temperatura inicial):
    T_centro(t) - T0 = Qmax * g(t)  =>  Qmax_limite = Qref * (T_LIM - T0) / max_t[Tc_ref(t) - T0].
Se verifica con una nueva corrida. El resultado es una APROXIMACION NUMERICA (Euler explicito
en el tiempo, volumenes finitos en el espacio): no es una garantia exacta del limite continuo.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import csv
import math
import sys

# Carpeta del ejercicio: todas las rutas parten de aqui, sin importar el directorio actual.
CARPETA = Path(__file__).resolve().parent
CARPETA_RESULTADOS = CARPETA / "resultados"

# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES (unidades: m, s, °C, W/(m K), kg/m³, J/(kg K), W/m³)
# ---------------------------------------------------------------------------
# Visualizacion
MOSTRAR_GRAFICOS = True     # ventana con Q(t) y T_centro(t)
MOSTRAR_ANIMACION = True    # ventana animada (mapa + T del centro vs t)
GUARDAR_ANIMACION = True    # resultados/evolucion_fuente_variable_2d.gif (Pillow)
NUM_FOTOGRAMAS = 100        # fotogramas equiespaciados en [0, TIEMPO_FINAL_S] (>= 2)
FPS_ANIMACION = 15          # SOLO velocidad de reproduccion; no cambia dt ni los tiempos fisicos
EJECUTAR_SENSIBILIDAD = True  # refinamiento de malla y de dt (tarda unos segundos mas)
# Sin ventanas (MOSTRAR_GRAFICOS y MOSTRAR_ANIMACION False) se usa el backend Agg.

# Geometria y malla
LADO_X_M, LADO_Y_M = 4.0, 2.0
CELDAS_X, CELDAS_Y = 80, 40            # pares (la sensibilidad usa la mitad y el doble)
PUNTO_CENTRO_M = (LADO_X_M / 2, LADO_Y_M / 2)   # punto fisico vigilado [m]
PROFUNDIDAD_M = 1.0                    # espesor fuera del plano [m] (energias en J por esa profundidad)

# Material homogeneo (mismos valores que scripts/demo_difusion_2d.py)
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0

# Temperaturas
TEMPERATURA_INICIAL_C = 20.0
TEMPERATURA_BORDE_C = 20.0             # Dirichlet en las cuatro caras (= T inicial: problema lineal)
TEMPERATURA_LIMITE_C = 70.0            # el centro no debe superarla

# Tiempos (tres cosas distintas)
DURACION_FUENTE_S = 2 * 24 * 3600.0    # tau: duracion de la fuente = 172800 s
TIEMPO_FINAL_S = 48 * 3600.0           # horizonte de evaluacion / fin de la simulacion
FRACCION_PASO_ESTABLE = 0.8            # el dt nunca supera esta fraccion del limite estable

QMAX_REFERENCIA_W_M3 = 100.0           # fuente de referencia positiva para escalar [W/m³]
# ---------------------------------------------------------------------------

import numpy as np  # noqa: E402

import matplotlib  # noqa: E402
if not (MOSTRAR_GRAFICOS or MOSTRAR_ANIMACION):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402


def cargar_solver():
    """Importa EXPLICITAMENTE solver_ejercicio.py de esta carpeta (no el del repositorio)."""
    ruta = CARPETA / "solver_ejercicio.py"
    spec = spec_from_file_location("solver_ejercicio", ruta)
    modulo = module_from_spec(spec)
    sys.modules["solver_ejercicio"] = modulo
    spec.loader.exec_module(modulo)
    return modulo


solver = cargar_solver()
_ANIMACIONES = []   # referencia viva de la animacion


# ---------------------------------------------------------------------------
# Fuente y utilidades
# ---------------------------------------------------------------------------
def fuente_q(t_s, qmax_w_m3, tau_s=None):
    """Q(t) = Qmax/2 [1 - cos(2 pi t/tau)] para 0 <= t <= tau; 0 fuera [W/m³]."""
    tau = DURACION_FUENTE_S if tau_s is None else tau_s
    if 0.0 <= t_s <= tau:
        return 0.5 * qmax_w_m3 * (1.0 - math.cos(2.0 * math.pi * t_s / tau))
    return 0.0


def energia_fuente_teorica_j(qmax_w_m3, volumen_total_m3, horizonte_s, tau_s=None):
    """Integral exacta de Q(t) V sobre [0, min(horizonte, tau)] [J]."""
    tau = DURACION_FUENTE_S if tau_s is None else tau_s
    h = min(horizonte_s, tau)
    return qmax_w_m3 * volumen_total_m3 * 0.5 * (h - tau / (2 * math.pi) * math.sin(2 * math.pi * h / tau))


def qmax_por_linealidad(qref_w_m3, exceso_max_ref_c, limite_c=None, t0_c=None):
    """Qmax_limite = Qref (T_lim - T0) / max_t[Tc_ref(t) - T0]."""
    limite = TEMPERATURA_LIMITE_C if limite_c is None else limite_c
    t0 = TEMPERATURA_INICIAL_C if t0_c is None else t0_c
    if exceso_max_ref_c <= 0:
        raise ValueError("la fuente de referencia no calienta el centro: no se puede escalar")
    return qref_w_m3 * (limite - t0) / exceso_max_ref_c


def construir_sistema(celdas_x, celdas_y, qmax_w_m3):
    fuente = lambda t: fuente_q(t, qmax_w_m3)  # noqa: E731
    bordes = {c: solver.dirichlet(TEMPERATURA_BORDE_C) for c in ("x-", "x+", "y-", "y+")}
    return solver.DifusionCartesiana(
        forma=(celdas_x, celdas_y), dx_m=(LADO_X_M / celdas_x, LADO_Y_M / celdas_y),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K, densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K, fuente_w_m3=fuente,
        fronteras=bordes, profundidad_m=PROFUNDIDAD_M)


def correr(celdas_x, celdas_y, qmax_w_m3, dt_s, capturar=False):
    """Una integracion continua; registra el centro en TODOS los pasos."""
    sistema = construir_sistema(celdas_x, celdas_y, qmax_w_m3)
    instantes = None
    if capturar:
        instantes = [j * TIEMPO_FINAL_S / (NUM_FOTOGRAMAS - 1) for j in range(NUM_FOTOGRAMAS)]
    resultado = sistema.simular(
        temperatura_inicial_c=TEMPERATURA_INICIAL_C, tiempo_final_s=TIEMPO_FINAL_S,
        paso_tiempo_max_s=dt_s * (1.0 + 1e-9), instantes_captura_s=instantes,
        puntos_registro_m=[PUNTO_CENTRO_M])
    return sistema, resultado


def maximo_del_centro(resultado):
    """(T maxima [°C], instante [s], indice de paso) sobre TODOS los pasos."""
    tc = resultado.historial_puntos_c[:, 0]
    k = int(np.argmax(tc))
    return float(tc[k]), float(resultado.historial_tiempos_s[k]), k


def paso_de_animacion():
    """dt tal que cada fotograma contiene un numero entero de pasos y dt <= fraccion*limite."""
    if NUM_FOTOGRAMAS < 2:
        raise ValueError("NUM_FOTOGRAMAS debe ser >= 2")
    sistema = construir_sistema(CELDAS_X, CELDAS_Y, 1.0)
    intervalo = TIEMPO_FINAL_S / (NUM_FOTOGRAMAS - 1)
    objetivo = FRACCION_PASO_ESTABLE * sistema.paso_maximo_estable_s()
    pasos = math.ceil(intervalo / objetivo)
    return intervalo / pasos, pasos


def sensibilidades(dt_base_s, qref):
    """Qmax_limite y t_max por separado: (a) malla con dt comun, (b) dt con malla fija."""
    filas_malla, filas_dt = [], []
    fino = construir_sistema(2 * CELDAS_X, 2 * CELDAS_Y, 1.0)
    k = math.ceil(dt_base_s / (FRACCION_PASO_ESTABLE * fino.paso_maximo_estable_s()))
    dt_comun = dt_base_s / k                      # estable en las tres mallas
    for factor in (0.5, 1, 2):
        cx, cy = int(CELDAS_X * factor), int(CELDAS_Y * factor)
        _, r = correr(cx, cy, qref, dt_comun)
        tmax, tiempo, _ = maximo_del_centro(r)
        filas_malla.append((f"{cx}x{cy}", r.paso_tiempo_s, qmax_por_linealidad(qref, tmax - TEMPERATURA_INICIAL_C), tmax, tiempo))
    for divisor in (1, 2, 4):
        _, r = correr(CELDAS_X, CELDAS_Y, qref, dt_base_s / divisor)
        tmax, tiempo, _ = maximo_del_centro(r)
        filas_dt.append((f"dt/{divisor}", r.paso_tiempo_s, qmax_por_linealidad(qref, tmax - TEMPERATURA_INICIAL_C), tmax, tiempo))
    return dt_comun, filas_malla, filas_dt


def tabla_sensibilidad(titulo, filas, etiqueta):
    ref = filas[-1]       # el mas fino / el dt mas pequeno
    L = [titulo, f"  {etiqueta:<10}{'dt [s]':>10}{'Qmax_lim [W/m³]':>18}{'dQmax vs ref':>14}{'t_max [h]':>11}{'dt_max vs ref [h]':>19}"]
    for nombre, dt, q, tmax, tiempo in filas:
        L.append(f"  {nombre:<10}{dt:>10.2f}{q:>18.5f}{(q - ref[2]) / ref[2]:>+14.3%}{tiempo / 3600:>11.3f}{(tiempo - ref[4]) / 3600:>+19.3f}")
    L.append("  (referencia de la variacion: ultima fila)")
    return L


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------
def guardar_csv(sistema_v, r_v, qmax):
    ruta = CARPETA_RESULTADOS / "centro_cada_paso.csv"
    t, tc = r_v.historial_tiempos_s, r_v.historial_puntos_c[:, 0]
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["paso", "tiempo_s", "tiempo_h", "Q_inicio_de_paso_W_m3", "T_centro_C"])
        for n, (ti, ci) in enumerate(zip(t, tc)):
            w.writerow([n, repr(float(ti)), repr(float(ti) / 3600.0), repr(fuente_q(float(ti), qmax)), repr(float(ci))])
    return ruta


def grafico_fuente(qmax, r_v, tmax, tiempo_max):
    fig, (ax_q, ax_t) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    tt = np.linspace(0.0, max(DURACION_FUENTE_S, TIEMPO_FINAL_S) * 1.1, 600)
    ax_q.plot(tt / 3600, [fuente_q(float(x), qmax) for x in tt], color="tab:red")
    ax_q.axvline(DURACION_FUENTE_S / 3600, color="gray", ls="--", label=f"tau = {DURACION_FUENTE_S / 3600:g} h")
    ax_q.axvline(TIEMPO_FINAL_S / 3600, color="k", ls=":", label=f"horizonte = {TIEMPO_FINAL_S / 3600:g} h")
    ax_q.set_ylabel("Q(t) [W/m³]"); ax_q.set_title(f"Fuente Q(t), Qmax = {qmax:.4f} W/m³"); ax_q.legend(); ax_q.grid(alpha=0.3)
    ax_t.plot(r_v.historial_tiempos_s / 3600, r_v.historial_puntos_c[:, 0], color="tab:blue", label="T centro")
    ax_t.axhline(TEMPERATURA_LIMITE_C, color="tab:red", ls="--", label=f"limite {TEMPERATURA_LIMITE_C:g} °C")
    ax_t.plot([tiempo_max / 3600], [tmax], "o", color="k", label=f"maximo {tmax:.3f} °C en {tiempo_max / 3600:.2f} h")
    ax_t.set_xlabel("Tiempo [h]"); ax_t.set_ylabel("T en el centro [°C]"); ax_t.legend(); ax_t.grid(alpha=0.3)
    fig.savefig(CARPETA_RESULTADOS / "grafico_fuente_Q_t.png", dpi=110)
    return fig


def crear_animacion(r_v, qmax, dt):
    campos, instantes = r_v.estados_capturados_c, np.array(r_v.instantes_capturados_s)
    pasos_frame = [int(round(t / dt)) for t in instantes]     # indice de paso de cada fotograma
    horas_hist = r_v.historial_tiempos_s / 3600.0
    tc = r_v.historial_puntos_c[:, 0]
    vmin = min(TEMPERATURA_BORDE_C, float(min(c.min() for c in campos)))
    vmax = float(max(TEMPERATURA_LIMITE_C, max(c.max() for c in campos)))
    fig, (eje, eje_t) = plt.subplots(1, 2, figsize=(13.0, 5.4), constrained_layout=True,
                                     gridspec_kw={"width_ratios": [1.25, 1.0]})
    imagen = eje.imshow(campos[0].T, origin="lower", extent=(0, LADO_X_M, 0, LADO_Y_M),
                        vmin=vmin, vmax=vmax, cmap="inferno", aspect="equal")
    fig.colorbar(imagen, ax=eje, label="T [°C]", fraction=0.046, pad=0.04)
    eje.plot(*PUNTO_CENTRO_M, "o", ms=9, mfc="tab:cyan", mec="white", mew=1.5)
    eje.set_xlabel("x [m]"); eje.set_ylabel("y [m]")
    eje_t.set_xlim(0.0, horas_hist[-1]); eje_t.set_ylim(vmin - 1, max(vmax, TEMPERATURA_LIMITE_C) + 2)
    eje_t.axhline(TEMPERATURA_LIMITE_C, color="tab:red", ls="--", lw=1, label=f"limite {TEMPERATURA_LIMITE_C:g} °C")
    linea, = eje_t.plot([], [], color="tab:cyan", lw=2, label="T centro")
    marca, = eje_t.plot([], [], "o", color="tab:cyan", ms=8, mec="white")
    eje_t.set_xlabel("Tiempo [h]"); eje_t.set_ylabel("Temperatura del centro [°C]"); eje_t.grid(alpha=0.3); eje_t.legend(loc="upper left")
    tit_mapa, tit_t = eje.set_title(""), eje_t.set_title("")

    def actualizar(i):
        k = pasos_frame[i]
        imagen.set_data(campos[i].T)
        texto = f"t = {instantes[i] / 3600:.2f} h, Q = {fuente_q(float(instantes[i]), qmax):.4f} W/m³"
        tit_mapa.set_text("Mapa de temperatura: " + texto)
        tit_t.set_text("Centro: " + texto)
        linea.set_data(horas_hist[:k + 1], tc[:k + 1])      # todos los pasos hasta este fotograma
        marca.set_data([horas_hist[k]], [tc[k]])
        return imagen, tit_mapa, tit_t, linea, marca

    actualizar(0)
    anim = FuncAnimation(fig, actualizar, frames=len(campos), interval=1000.0 / FPS_ANIMACION, blit=False, repeat=True)
    _ANIMACIONES.append(anim)
    if GUARDAR_ANIMACION:
        ruta = CARPETA_RESULTADOS / "evolucion_fuente_variable_2d.gif"
        anim.save(ruta, writer=PillowWriter(fps=FPS_ANIMACION), dpi=80)
        print(f"Animacion guardada: {ruta}")
    if not MOSTRAR_ANIMACION:
        plt.close(fig)
    return anim


def main():
    CARPETA_RESULTADOS.mkdir(parents=True, exist_ok=True)
    qref = QMAX_REFERENCIA_W_M3
    dt_base, pasos_frame = paso_de_animacion()

    # 1) Referencia positiva -> Qmax_limite por linealidad (maximo sobre TODOS los pasos).
    _, r_ref = correr(CELDAS_X, CELDAS_Y, qref, dt_base)
    tmax_ref, t_ref, _ = maximo_del_centro(r_ref)
    qmax = qmax_por_linealidad(qref, tmax_ref - TEMPERATURA_INICIAL_C)

    # 2) Verificacion con una NUEVA corrida con Qmax_limite (tambien produce la animacion).
    sistema_v, r_v = correr(CELDAS_X, CELDAS_Y, qmax, dt_base, capturar=True)
    tmax, tiempo_max, k_max = maximo_del_centro(r_v)
    tc = r_v.historial_puntos_c[:, 0]

    sens = sensibilidades(dt_base, qref) if EJECUTAR_SENSIBILIDAD else None

    # 3) Resumen.
    vol_total = sistema_v.volumen_m3 * CELDAS_X * CELDAS_Y
    cambio = r_v.energia_final_j - r_v.energia_inicial_j
    e_teo = energia_fuente_teorica_j(qmax, vol_total, TIEMPO_FINAL_S)
    L = ["EJERCICIO FUENTE VARIABLE 2D - RESUMEN", "",
         "RESULTADO (aproximacion numerica, no garantia exacta del limite continuo)",
         f"  Qmax calculado (limite): {qmax:.6f} W/m³  (referencia {qref:g} W/m³ -> centro max {tmax_ref:.4f} °C)",
         f"  Temperatura maxima del centro (corrida de verificacion): {tmax:.6f} °C  (limite {TEMPERATURA_LIMITE_C:g} °C, diferencia {tmax - TEMPERATURA_LIMITE_C:+.2e} °C)",
         f"  Instante del maximo: {tiempo_max:.1f} s = {tiempo_max / 3600:.3f} h (paso {k_max} de {r_v.numero_pasos}; horizonte {TIEMPO_FINAL_S / 3600:g} h)",
         f"  Temperatura del centro {PUNTO_CENTRO_M} m: inicial {tc[0]:.3f} °C, final {tc[-1]:.3f} °C",
         "  El maximo se busca en TODOS los pasos del solver, no solo en fotogramas ni al final.",
         "", "CASO FISICO",
         f"  Rectangulo {LADO_X_M:g} m x {LADO_Y_M:g} m (profundidad {PROFUNDIDAD_M:g} m); T inicial {TEMPERATURA_INICIAL_C:g} °C; cuatro bordes Dirichlet a {TEMPERATURA_BORDE_C:g} °C",
         f"  Propiedades (homogeneas): k = {CONDUCTIVIDAD_W_M_K:g} W/(m K), rho = {DENSIDAD_KG_M3:g} kg/m³, cp = {CALOR_ESPECIFICO_J_KG_K:g} J/(kg K), alfa = {CONDUCTIVIDAD_W_M_K / (DENSIDAD_KG_M3 * CALOR_ESPECIFICO_J_KG_K):.3e} m²/s",
         f"  Fuente: Q(t) = Qmax/2 [1 - cos(2 pi t/tau)], tau = {DURACION_FUENTE_S:g} s = {DURACION_FUENTE_S / 3600:g} h; tiempo final = {TIEMPO_FINAL_S:g} s = {TIEMPO_FINAL_S / 3600:g} h; FPS solo reproduccion = {FPS_ANIMACION}",
         "", "MALLA Y PASO",
         f"  Celdas {CELDAS_X} x {CELDAS_Y}; dx = {LADO_X_M / CELDAS_X:g} m, dy = {LADO_Y_M / CELDAS_Y:g} m",
         f"  dt efectivo = {r_v.paso_tiempo_s:.3f} s; limite estable = {r_v.paso_maximo_estable_s:.3f} s ({r_v.paso_tiempo_s / r_v.paso_maximo_estable_s:.1%}); pasos = {r_v.numero_pasos} ({pasos_frame} por fotograma, {NUM_FOTOGRAMAS} fotogramas)",
         "  La fuente se evalua al inicio de cada paso: Q(t_n) con t_n = n dt (Euler explicito).",
         "", "BALANCE DE ENERGIA (referencia 0 °C; frontera positiva = entra al cuerpo, negativa = sale; J por la profundidad usada)",
         f"  Energia inicial {r_v.energia_inicial_j:.6e} J; final {r_v.energia_final_j:.6e} J; cambio {cambio:+.6e} J",
         f"  Aporte de la fuente (aplicado paso a paso): {r_v.energia_fuente_j:+.6e} J  (integral exacta de Q: {e_teo:+.6e} J; dif. {r_v.energia_fuente_j - e_teo:+.2e} J = error de integracion temporal)",
         f"  Intercambio neto de fronteras: {r_v.energia_frontera_total_j:+.6e} J",
         f"  Residuo = cambio - fuente - fronteras: {r_v.error_balance_j:+.3e} J"]
    if sens:
        dt_comun, f_malla, f_dt = sens
        L += ["", "SENSIBILIDAD (cada una por separado; Qmax_limite por linealidad con Qref)"]
        L += tabla_sensibilidad(f"  (a) Refinamiento de malla con dt comun estable = {dt_comun:.2f} s", f_malla, "malla")
        L += tabla_sensibilidad(f"  (b) Reduccion de dt con malla fija {CELDAS_X}x{CELDAS_Y}", f_dt, "paso")
        L += ["  Las variaciones son indicadores de convergencia numerica, no cotas de error."]
    texto = "\n".join(L)
    (CARPETA_RESULTADOS / "resumen.txt").write_text(texto + "\n", encoding="utf-8")
    print(texto)
    print(f"Resumen guardado: {CARPETA_RESULTADOS / 'resumen.txt'}")
    print(f"CSV guardado: {guardar_csv(sistema_v, r_v, qmax)}")

    fig = grafico_fuente(qmax, r_v, tmax, tiempo_max)
    print(f"Grafico guardado: {CARPETA_RESULTADOS / 'grafico_fuente_Q_t.png'}")
    if not MOSTRAR_GRAFICOS:
        plt.close(fig)
    if MOSTRAR_ANIMACION or GUARDAR_ANIMACION:
        crear_animacion(r_v, qmax, r_v.paso_tiempo_s)
    if MOSTRAR_GRAFICOS or MOSTRAR_ANIMACION:
        plt.show()      # una sola llamada para todas las ventanas


if __name__ == "__main__":
    main()
