r"""Ejercicio de clases 3D: caja con fuente uniforme variable en el tiempo y distintas condiciones de borde.

Pulsa "Ejecutar archivo de Python" (VS Code) o:  python ejercicio_fuente_variable_3d.py
Funciona desde cualquier directorio. Experimento LOCAL, separado de E04/E05 y de los ejercicios 2D:
usa su propia copia del solver (solver_ejercicio.py, en esta carpeta).

Caso: caja Lx x Ly x Lz = 3 x 3 x 0.60 m, T inicial uniforme 20 °C en TODOS los escenarios, material homogeneo
y fuente uniforme en el espacio
    Q(t) = Qmax/2 [1 - cos(2 pi t / tau)]   (0 <= t <= tau),    Q(t) = 0 fuera.
Qmax NO se recalcula: es el valor del ejercicio Dirichlet 2D. Escenarios (cada uno es una simulacion aparte):
  1. Dirichlet 20 °C en las seis caras.
  2. Dirichlet 30 °C en las seis caras (interior inicial a 20 °C; el borde vale 30 °C desde t = 0).
  3. Neumann con flujo cero (aislamiento).
  4. Robin en las seis caras, ambiente a 20 °C, para varios beta.

CONVENCION ROBIN. Con beta = h/k > 0 [1/m] y n_ext la normal EXTERIOR:
    -k dT/dn_ext = h (T_superficie - T_ambiente)   <=>   dT/dn_ext = -beta (T_superficie - T_ambiente).
Si el cuerpo esta mas caliente que el ambiente, el flujo sale (se enfria); si esta mas frio, entra. El pizarron
muestra un signo positivo: eso corresponde a derivar respecto de la normal INTERIOR (dT/dn_int = +beta (T_s - T_amb)),
que es lo mismo. El solver usa robin(h, T_ambiente), por eso h = k * beta. Los beta de prueba son EXPLORATORIOS,
no datos entregados por el profesor.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import csv
import math
import sys
import time

CARPETA = Path(__file__).resolve().parent
CARPETA_RESULTADOS = CARPETA / "resultados"

# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES (unidades: m, s, °C, W/(m K), kg/m³, J/(kg K), W/m³, 1/m)
# ---------------------------------------------------------------------------
# Visualizacion
MOSTRAR_GRAFICOS = True     # ventana con la comparacion de escenarios (Q(t) y T central)
MOSTRAR_ANIMACION = True    # ventanas animadas (una por escenario animado)
GUARDAR_ANIMACION = True    # un GIF por escenario en resultados/<escenario>/
ESCENARIOS_ANIMADOS = None  # None = todos; o una lista de nombres, p. ej. ["neumann_cero", "robin_beta_1"]
NUM_FOTOGRAMAS = 60         # fotogramas equiespaciados en [0, TIEMPO_FINAL_S] (>= 2)
FPS_ANIMACION = 12          # SOLO velocidad de reproduccion; no cambia dt ni los tiempos fisicos
EJECUTAR_SENSIBILIDAD = True  # refinamiento de malla y de dt (mas lento; ver el costo en el resumen)
# Sin ventanas (MOSTRAR_GRAFICOS y MOSTRAR_ANIMACION False) se usa el backend Agg.

# Geometria y malla (dx = dy = dz = 0.05 m)
LADO_X_M, LADO_Y_M, LADO_Z_M = 3.0, 3.0, 0.60
CELDAS_X, CELDAS_Y, CELDAS_Z = 60, 60, 12
PUNTO_CENTRO_M = (LADO_X_M / 2, LADO_Y_M / 2, LADO_Z_M / 2)     # punto fisico vigilado [m]

# Material homogeneo: mismos valores que el ejercicio de fuente variable 2D
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0

# Fuente. Qmax FIJO (no se recalcula ni se optimiza para esta geometria).
# Procedencia: ejercicio Dirichlet 2D (ejercicio_fuente_variable_2d.py, malla 80x40, dt = 349.09 s). El resumen
# de ese ejercicio publica 1282.565913; el valor completo se recupero de
# ejercicios_locales/fuente_variable_2d/resultados/centro_cada_paso.csv despejando Q(t_n) = Qmax/2 [1 - cos(2 pi t_n/tau)]
# en 312 pasos (dispersion 5e-14 W/m³).
QMAX_W_M3 = 1282.5659131507364
DURACION_FUENTE_S = 2 * 24 * 3600.0    # tau = 172800 s (independiente del horizonte)
TIEMPO_FINAL_S = 48 * 3600.0           # horizonte de simulacion (independiente de tau)

# Temperaturas
TEMPERATURA_INICIAL_C = 20.0           # interior, en todos los escenarios
TEMPERATURA_DIRICHLET_1_C = 20.0       # escenario 1
TEMPERATURA_DIRICHLET_2_C = 30.0       # escenario 2
TEMPERATURA_AMBIENTE_C = 20.0          # escenario 4 (Robin)
TEMPERATURA_REFERENCIA_C = 70.0        # linea de referencia: superarla es un resultado a informar, no a corregir

# Robin: beta = h/k [1/m]. VALORES EXPLORATORIOS (no son datos del profesor).
BETAS_ROBIN_1_M = [0.1, 1.0, 10.0]

FRACCION_PASO_ESTABLE = 0.8            # el dt comun nunca supera esta fraccion del limite del escenario mas exigente
# Sensibilidad: mallas = factor * (CELDAS_X, CELDAS_Y, CELDAS_Z); dt = dt_base / divisor
SENS_ESCENARIOS = ["dirichlet_20C", "robin_beta_1"]
SENS_FACTORES_MALLA = (0.5, 1, 2)
SENS_DIVISORES_DT = (1, 2, 4)
# ---------------------------------------------------------------------------

import numpy as np  # noqa: E402

import matplotlib  # noqa: E402
if not (MOSTRAR_GRAFICOS or MOSTRAR_ANIMACION):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402


def cargar_solver():
    """Importa EXPLICITAMENTE solver_ejercicio.py de esta carpeta (no el del repositorio ni el del 2D)."""
    spec = spec_from_file_location("solver_ejercicio_3d", CARPETA / "solver_ejercicio.py")
    modulo = module_from_spec(spec)
    sys.modules["solver_ejercicio_3d"] = modulo
    spec.loader.exec_module(modulo)
    return modulo


solver = cargar_solver()
CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")
RHOCP = DENSIDAD_KG_M3 * CALOR_ESPECIFICO_J_KG_K
_ANIMACIONES = []


# ---------------------------------------------------------------------------
# Fuente, escenarios y solucion adiabatica exacta
# ---------------------------------------------------------------------------
def fuente_q(t_s, qmax_w_m3=None):
    """Q(t) = Qmax/2 [1 - cos(2 pi t/tau)] para 0 <= t <= tau; 0 fuera [W/m³]."""
    qmax = QMAX_W_M3 if qmax_w_m3 is None else qmax_w_m3
    if 0.0 <= t_s <= DURACION_FUENTE_S:
        return 0.5 * qmax * (1.0 - math.cos(2.0 * math.pi * t_s / DURACION_FUENTE_S))
    return 0.0


def h_desde_beta(beta_1_m, k_w_m_k=None):
    """h = k beta [W/(m² K)] con beta = h/k [1/m]."""
    return (CONDUCTIVIDAD_W_M_K if k_w_m_k is None else k_w_m_k) * beta_1_m


def nombre_beta(beta):
    return "robin_beta_" + f"{beta:g}".replace(".", "p")


def definir_escenarios():
    """Lista de escenarios: dict(slug, titulo, fronteras, beta, h)."""
    esc = [
        dict(slug="dirichlet_20C", titulo=f"Dirichlet {TEMPERATURA_DIRICHLET_1_C:g} °C", beta=None, h=None,
             fronteras={c: solver.dirichlet(TEMPERATURA_DIRICHLET_1_C) for c in CARAS}),
        dict(slug="dirichlet_30C", titulo=f"Dirichlet {TEMPERATURA_DIRICHLET_2_C:g} °C", beta=None, h=None,
             fronteras={c: solver.dirichlet(TEMPERATURA_DIRICHLET_2_C) for c in CARAS}),
        dict(slug="neumann_cero", titulo="Neumann q = 0 (aislado)", beta=None, h=None,
             fronteras={c: solver.aislado() for c in CARAS}),
    ]
    for beta in BETAS_ROBIN_1_M:
        h = h_desde_beta(beta)
        esc.append(dict(slug=nombre_beta(beta), titulo=f"Robin beta = {beta:g} 1/m (h = {h:g} W/(m² K))",
                        beta=beta, h=h, fronteras={c: solver.robin(h, TEMPERATURA_AMBIENTE_C) for c in CARAS}))
    return esc


def temperatura_adiabatica_exacta(t_s, qmax_w_m3=None):
    """T(t) de un dominio aislado con fuente uniforme y T inicial uniforme [°C] (vectorizada)."""
    qmax = QMAX_W_M3 if qmax_w_m3 is None else qmax_w_m3
    tau = DURACION_FUENTE_S
    t = np.minimum(np.asarray(t_s, dtype=float), tau)
    return TEMPERATURA_INICIAL_C + qmax / (2.0 * RHOCP) * (t - tau / (2.0 * math.pi) * np.sin(2.0 * math.pi * t / tau))


def tolerancia_integracion_c(dt_s, qmax_w_m3=None):
    """Cota del error de integracion temporal (suma por la izquierda, Euler explicito); ver el ejercicio 2D:
    (dt/2) Qmax/(rho cp) + (dt²/12)(pi Qmax/tau)/(rho cp)."""
    qmax = QMAX_W_M3 if qmax_w_m3 is None else qmax_w_m3
    return (0.5 * dt_s * qmax + dt_s**2 / 12.0 * math.pi * qmax / DURACION_FUENTE_S) / RHOCP


# ---------------------------------------------------------------------------
# Simulacion
# ---------------------------------------------------------------------------
def construir_sistema(fronteras, celdas=None, fuente=None, lados=None):
    celdas = celdas or (CELDAS_X, CELDAS_Y, CELDAS_Z)
    lados = lados or (LADO_X_M, LADO_Y_M, LADO_Z_M)
    return solver.DifusionCartesiana(
        forma=celdas, dx_m=tuple(L / n for L, n in zip(lados, celdas)),
        conductividad_w_m_k=CONDUCTIVIDAD_W_M_K, densidad_kg_m3=DENSIDAD_KG_M3,
        calor_especifico_j_kg_k=CALOR_ESPECIFICO_J_KG_K,
        fuente_w_m3=(lambda t: fuente_q(t)) if fuente is None else fuente, fronteras=fronteras)


def paso_comun(escenarios):
    """dt COMUN estable para todos los escenarios (caras, aristas y esquinas incluidas en el limite del solver).

    dt <= FRACCION * min(limites); ademas cada fotograma contiene un numero entero m de pasos."""
    limites = {e["slug"]: construir_sistema(e["fronteras"]).paso_maximo_estable_s() for e in escenarios}
    intervalo = TIEMPO_FINAL_S / (NUM_FOTOGRAMAS - 1)
    objetivo = FRACCION_PASO_ESTABLE * min(limites.values())
    pasos = math.ceil(intervalo / objetivo)
    return intervalo / pasos, pasos, limites


def simular(esc, dt_s, celdas=None, capturar=False, fuente=None, T0=None, tiempo_final_s=None):
    sistema = construir_sistema(esc["fronteras"], celdas, fuente)
    tf = TIEMPO_FINAL_S if tiempo_final_s is None else tiempo_final_s
    instantes = [j * tf / (NUM_FOTOGRAMAS - 1) for j in range(NUM_FOTOGRAMAS)] if capturar else None
    t0 = time.perf_counter()
    r = sistema.simular(temperatura_inicial_c=TEMPERATURA_INICIAL_C if T0 is None else T0, tiempo_final_s=tf,
                        paso_tiempo_max_s=dt_s * (1.0 + 1e-9), instantes_captura_s=instantes,
                        puntos_registro_m=[PUNTO_CENTRO_M], registrar_metricas=True)
    return sistema, r, time.perf_counter() - t0


def corte_interpolado(T, eje, coord_m, dx_m):
    """Plano T(eje = coord) por interpolacion LINEAL entre capas de centros. Devuelve (plano, descripcion)."""
    n, d = T.shape[eje], dx_m[eje]
    s = min(max(coord_m / d - 0.5, 0.0), n - 1.0)
    i = min(int(math.floor(s)), max(n - 2, 0))
    w = s - i
    plano = (1.0 - w) * np.take(T, i, axis=eje) + w * np.take(T, min(i + 1, n - 1), axis=eje)
    c0, c1 = (i + 0.5) * d, (min(i + 1, n - 1) + 0.5) * d
    desc = (f"{'xyz'[eje]} = {coord_m:.3f} m (capa {c0:.3f} m)" if w == 0.0 else
            f"{'xyz'[eje]} = {coord_m:.3f} m (interpolado entre capas {c0:.3f} y {c1:.3f} m, pesos {1 - w:.2f}/{w:.2f})")
    return plano, desc


def tres_cortes(T, dx_m):
    xy, d_xy = corte_interpolado(T, 2, PUNTO_CENTRO_M[2], dx_m)
    xz, d_xz = corte_interpolado(T, 1, PUNTO_CENTRO_M[1], dx_m)
    yz, d_yz = corte_interpolado(T, 0, PUNTO_CENTRO_M[0], dx_m)
    return (xy, xz, yz), (d_xy, d_xz, d_yz)


def primer_cruce_s(tiempos, valores, limite):
    k = np.flatnonzero(valores >= limite)
    if k.size == 0:
        return None
    k = int(k[0])
    if k == 0:
        return float(tiempos[0])
    t0, t1, v0, v1 = tiempos[k - 1], tiempos[k], valores[k - 1], valores[k]
    return float(t0 + (limite - v0) * (t1 - t0) / (v1 - v0))


def simetria_max_dif(T):
    """Mayor diferencia entre T y sus reflexiones en x, y, z y el intercambio x <-> y (Lx = Ly) [°C]."""
    return max(float(np.max(np.abs(T - T[::-1, :, :]))), float(np.max(np.abs(T - T[:, ::-1, :]))),
               float(np.max(np.abs(T - T[:, :, ::-1]))), float(np.max(np.abs(T - T.transpose(1, 0, 2)))))


def metricas(esc, sistema, r, dt, segundos):
    t, tc = r.historial_tiempos_s, r.historial_puntos_c[:, 0]
    k_c = int(np.argmax(tc)); k_g = int(np.argmax(r.historial_max_c))
    return dict(
        slug=esc["slug"], titulo=esc["titulo"], beta=esc["beta"], h=esc["h"],
        Tc_ini=float(tc[0]), Tc_fin=float(tc[-1]), Tc_max=float(tc[k_c]), t_Tc_max_s=float(t[k_c]),
        Tg_max=float(r.historial_max_c[k_g]), t_Tg_max_s=float(t[k_g]), Tg_min=float(r.historial_min_c.min()),
        cruce_s=primer_cruce_s(t, tc, TEMPERATURA_REFERENCIA_C), supera=bool(tc.max() > TEMPERATURA_REFERENCIA_C),
        dt=r.paso_tiempo_s, limite=r.paso_maximo_estable_s, pasos=r.numero_pasos, segundos=segundos,
        e_ini=r.energia_inicial_j, e_fin=r.energia_final_j, e_fuente=r.energia_fuente_j,
        e_front=r.energia_frontera_total_j, e_caras=dict(r.energia_frontera_j), residuo=r.error_balance_j,
        simetria=simetria_max_dif(r.temperatura_c), dispersion_final=float(r.temperatura_c.max() - r.temperatura_c.min()),
        celdas=sistema.forma)


# ---------------------------------------------------------------------------
# Verificaciones rapidas (malla gruesa) - las mismas funciones las usan los tests
# ---------------------------------------------------------------------------
CELDAS_VERIFICACION = (12, 12, 3)


def verificaciones_rapidas():
    """Devuelve lista de (nombre, aprobado, detalle) con casos de control en malla gruesa y 6 h."""
    out = []
    horizonte = 6 * 3600.0
    esc = {e["slug"]: e for e in definir_escenarios()}
    # dt estable para TODOS los casos de control (el Dirichlet es el mas restrictivo)
    dt = FRACCION_PASO_ESTABLE * min(construir_sistema(e["fronteras"], CELDAS_VERIFICACION).paso_maximo_estable_s()
                                     for e in esc.values())

    # a) Robin con h = 0 equivale a aislamiento
    robin0 = dict(slug="robin0", titulo="", beta=0.0, h=0.0, fronteras={c: solver.robin(0.0, TEMPERATURA_AMBIENTE_C) for c in CARAS})
    _, r_a, _ = simular(esc["neumann_cero"], dt, CELDAS_VERIFICACION, tiempo_final_s=horizonte)
    _, r_b, _ = simular(robin0, dt, CELDAS_VERIFICACION, tiempo_final_s=horizonte)
    dif = float(np.max(np.abs(r_a.temperatura_c - r_b.temperatura_c)))
    out.append(("Robin h = 0 equivale a aislamiento", dif <= 1e-12 and r_b.energia_frontera_total_j == 0.0, f"6 h; max|dT| = {dif:.2e} °C"))

    # b) Sin fuente, Robin enfria un cuerpo mas caliente que el ambiente (y no lo calienta)
    robin1 = dict(slug="r1", titulo="", beta=1.0, h=h_desde_beta(1.0), fronteras={c: solver.robin(h_desde_beta(1.0), TEMPERATURA_AMBIENTE_C) for c in CARAS})
    _, r_c, _ = simular(robin1, dt, CELDAS_VERIFICACION, fuente=lambda t: 0.0, T0=40.0, tiempo_final_s=horizonte)
    tc = r_c.historial_puntos_c[:, 0]
    ok = (r_c.energia_frontera_total_j < 0.0 and r_c.temperatura_c.max() < 40.0 and r_c.temperatura_c.min() >= TEMPERATURA_AMBIENTE_C
          and bool(np.all(np.diff(tc) <= 1e-12)))
    out.append(("Sin fuente, Robin enfria un cuerpo caliente", ok,
                f"6 h; T media final {r_c.temperatura_c.mean():.3f} °C (inicial 40), energia de frontera {r_c.energia_frontera_total_j:+.3e} J"))

    # c) Sin fuente y T inicial = Dirichlet: equilibrio
    _, r_d, _ = simular(esc["dirichlet_20C"], dt, CELDAS_VERIFICACION, fuente=lambda t: 0.0, tiempo_final_s=horizonte)
    dev = float(np.max(np.abs(r_d.temperatura_c - TEMPERATURA_DIRICHLET_1_C)))
    out.append(("Sin fuente y T inicial = Dirichlet: equilibrio", dev <= 1e-12, f"6 h; max|T - 20| = {dev:.2e} °C"))

    # d) Neumann cero: campo uniforme y solucion adiabatica analitica (horizonte completo)
    _, r_e, _ = simular(esc["neumann_cero"], dt / 8, CELDAS_VERIFICACION)      # dt pequeno: el error de Euler es chico
    err = float(np.max(np.abs(r_e.historial_puntos_c[:, 0] - temperatura_adiabatica_exacta(r_e.historial_tiempos_s))))
    tol = tolerancia_integracion_c(r_e.paso_tiempo_s)
    unif = float(r_e.temperatura_c.max() - r_e.temperatura_c.min())
    out.append(("Neumann cero: uniforme y solucion adiabatica", err <= tol and unif <= 1e-9 and r_e.energia_frontera_total_j == 0.0,
                f"{TIEMPO_FINAL_S / 3600:g} h, dt/8; error {err:.3e} <= tol {tol:.3e} °C; dispersion {unif:.2e} °C; T final {r_e.temperatura_c.mean():.4f} °C"))

    # e) Balance de energia con fuente y Robin
    _, r_f, _ = simular(robin1, dt, CELDAS_VERIFICACION)
    escala = abs(r_f.energia_fuente_j)
    out.append(("Balance de energia (Robin con fuente)", abs(r_f.error_balance_j) <= 1e-9 * escala, f"{TIEMPO_FINAL_S / 3600:g} h; residuo {r_f.error_balance_j:+.3e} J ({abs(r_f.error_balance_j) / escala:.1e} de la fuente)"))
    return out


# ---------------------------------------------------------------------------
# Salidas
# ---------------------------------------------------------------------------
def texto_escenario(m, dt_com, pasos_frame, cortes_desc):
    L = [f"ESCENARIO: {m['titulo']}  [{m['slug']}]", ""]
    if m["beta"] is not None:
        L.append(f"  Robin: beta = {m['beta']:g} 1/m -> h = k·beta = {m['h']:g} W/(m² K); ambiente {TEMPERATURA_AMBIENTE_C:g} °C; "
                 "-k dT/dn_ext = h (T_s - T_amb) (enfria si T_s > T_amb)")
    L += [f"  Caja {LADO_X_M:g} x {LADO_Y_M:g} x {LADO_Z_M:g} m, malla {m['celdas'][0]}x{m['celdas'][1]}x{m['celdas'][2]}; T inicial {TEMPERATURA_INICIAL_C:g} °C",
          f"  k = {CONDUCTIVIDAD_W_M_K:g} W/(m K), rho = {DENSIDAD_KG_M3:g} kg/m³, cp = {CALOR_ESPECIFICO_J_KG_K:g} J/(kg K), alfa = {CONDUCTIVIDAD_W_M_K / RHOCP:.3e} m²/s",
          f"  Fuente: Qmax = {QMAX_W_M3!r} W/m³, tau = {DURACION_FUENTE_S:g} s; horizonte = {TIEMPO_FINAL_S:g} s = {TIEMPO_FINAL_S / 3600:g} h",
          f"  dt efectivo (comun a todos) = {m['dt']:.3f} s; limite estable de este escenario = {m['limite']:.3f} s ({m['dt'] / m['limite']:.1%}); pasos = {m['pasos']} ({pasos_frame} por fotograma)",
          "", f"  Centro {PUNTO_CENTRO_M} m (interpolado): inicial {m['Tc_ini']:.4f} °C, final {m['Tc_fin']:.4f} °C, maximo {m['Tc_max']:.4f} °C en {m['t_Tc_max_s'] / 3600:.3f} h",
          f"  Maximo del DOMINIO (todas las celdas, todos los pasos): {m['Tg_max']:.4f} °C en {m['t_Tg_max_s'] / 3600:.3f} h (diferencia con el maximo central {m['Tg_max'] - m['Tc_max']:+.2e} °C); minimo del dominio {m['Tg_min']:.4f} °C",
          f"  Supera {TEMPERATURA_REFERENCIA_C:g} °C: {'SI' if m['supera'] else 'NO'}" + (f"; primer cruce (interpolacion lineal entre pasos) en {m['cruce_s'] / 3600:.3f} h" if m["cruce_s"] is not None else ""),
          f"  Simetria final (reflexiones x, y, z y x<->y): {m['simetria']:.2e} °C; dispersion espacial final max-min: {m['dispersion_final']:.3e} °C",
          "", "  Balance de energia (referencia 0 °C; frontera positiva = entra al cuerpo, negativa = sale; J)",
          f"    Inicial {m['e_ini']:.6e}; final {m['e_fin']:.6e}; cambio {m['e_fin'] - m['e_ini']:+.6e}",
          f"    Fuente aplicada {m['e_fuente']:+.6e}; fronteras netas {m['e_front']:+.6e} (" + ", ".join(f"{c} {v:+.3e}" for c, v in m["e_caras"].items()) + ")",
          f"    Residuo = cambio - fuente - fronteras: {m['residuo']:+.3e} J ({abs(m['residuo']) / max(abs(m['e_fuente']), 1e-300):.1e} de la fuente)",
          "", "  Cortes de la animacion: " + "; ".join(cortes_desc), f"  Tiempo de calculo: {m['segundos']:.1f} s"]
    return "\n".join(L)


def guardar_csv_escenario(carpeta, r):
    ruta = carpeta / "centro_cada_paso.csv"
    t = r.historial_tiempos_s
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["paso", "tiempo_s", "tiempo_h", "Q_inicio_de_paso_W_m3", "T_centro_C", "T_max_dominio_C", "T_min_dominio_C", "T_media_dominio_C"])
        for n in range(len(t)):
            w.writerow([n, repr(float(t[n])), repr(float(t[n]) / 3600.0), repr(fuente_q(float(t[n]))),
                        repr(float(r.historial_puntos_c[n, 0])), repr(float(r.historial_max_c[n])),
                        repr(float(r.historial_min_c[n])), repr(float(r.historial_media_c[n]))])
    return ruta


def guardar_tabla_comparativa(ms):
    ruta = CARPETA_RESULTADOS / "comparacion" / "tabla_comparativa.csv"
    ruta.parent.mkdir(parents=True, exist_ok=True)
    cols = ["escenario", "beta_1_m", "h_W_m2K", "Tc_inicial_C", "Tc_final_C", "Tc_max_C", "t_Tc_max_h", "Tglobal_max_C",
            "t_Tglobal_max_h", "supera_70C", "primer_cruce_70C_h", "dt_s", "residuo_balance_J"]
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for m in ms:
            w.writerow([m["slug"], "" if m["beta"] is None else m["beta"], "" if m["h"] is None else m["h"],
                        repr(m["Tc_ini"]), repr(m["Tc_fin"]), repr(m["Tc_max"]), repr(m["t_Tc_max_s"] / 3600),
                        repr(m["Tg_max"]), repr(m["t_Tg_max_s"] / 3600), m["supera"],
                        "" if m["cruce_s"] is None else repr(m["cruce_s"] / 3600), repr(m["dt"]), repr(m["residuo"])])
    return ruta


def tabla_texto(ms):
    L = [f"  {'escenario':<18}{'beta':>7}{'h':>8}{'Tc ini':>8}{'Tc fin':>9}{'Tc max':>9}{'t Tc max[h]':>12}{'Tglob max':>11}{'t Tglob[h]':>11}{'>70 C':>7}{'cruce70[h]':>11}"]
    for m in ms:
        L.append(f"  {m['slug']:<18}{'-' if m['beta'] is None else format(m['beta'], 'g'):>7}{'-' if m['h'] is None else format(m['h'], 'g'):>8}"
                 f"{m['Tc_ini']:>8.3f}{m['Tc_fin']:>9.3f}{m['Tc_max']:>9.3f}{m['t_Tc_max_s'] / 3600:>12.3f}{m['Tg_max']:>11.3f}"
                 f"{m['t_Tg_max_s'] / 3600:>11.3f}{'SI' if m['supera'] else 'no':>7}{'-' if m['cruce_s'] is None else format(m['cruce_s'] / 3600, '.3f'):>11}")
    L.append("  (beta en 1/m, h en W/(m² K), temperaturas en °C; Tc = centro, Tglob = maximo de todas las celdas)")
    return L


def sensibilidad(escenarios, dt_base):
    """Sensibilidad de malla (dt comun) y de dt (malla fija) para las metricas centrales de cada escenario elegido."""
    L = []
    por_slug = {e["slug"]: e for e in escenarios}
    for slug in SENS_ESCENARIOS:
        if slug not in por_slug:
            continue
        esc = por_slug[slug]
        mallas = [tuple(int(round(f * n)) for n in (CELDAS_X, CELDAS_Y, CELDAS_Z)) for f in SENS_FACTORES_MALLA]
        fino = construir_sistema(esc["fronteras"], mallas[-1])
        k = math.ceil(dt_base / (FRACCION_PASO_ESTABLE * fino.paso_maximo_estable_s()))
        dt_c = dt_base / k
        L += [f"  Escenario {slug}", f"    (a) Malla con dt comun estable = {dt_c:.2f} s",
              f"    {'malla':<12}{'Tc max':>10}{'dif vs fina':>15}{'t Tc max[h]':>13}{'Tc final':>10}{'Tglob max':>11}{'costo [s]':>11}"]
        filas = []
        for celdas in mallas:
            s, r, seg = simular(esc, dt_c, celdas)
            tc = r.historial_puntos_c[:, 0]; kk = int(np.argmax(tc))
            filas.append((celdas, float(tc[kk]), float(r.historial_tiempos_s[kk]), float(tc[-1]), float(r.historial_max_c.max()), seg))
        for celdas, tmax, tt, tf, tg, seg in filas:
            L.append(f"    {'x'.join(map(str, celdas)):<12}{tmax:>10.4f}{tmax - filas[-1][1]:>+15.4f}{tt / 3600:>13.3f}{tf:>10.4f}{tg:>11.4f}{seg:>11.1f}")
        L += [f"    (b) dt con malla fija {CELDAS_X}x{CELDAS_Y}x{CELDAS_Z}",
              f"    {'dt [s]':<12}{'Tc max':>10}{'dif vs dt min':>15}{'t Tc max[h]':>13}{'Tc final':>10}{'Tglob max':>11}{'costo [s]':>11}"]
        filas = []
        for div in SENS_DIVISORES_DT:
            s, r, seg = simular(esc, dt_base / div)
            tc = r.historial_puntos_c[:, 0]; kk = int(np.argmax(tc))
            filas.append((r.paso_tiempo_s, float(tc[kk]), float(r.historial_tiempos_s[kk]), float(tc[-1]), float(r.historial_max_c.max()), seg))
        for dtt, tmax, tt, tf, tg, seg in filas:
            L.append(f"    {dtt:<12.2f}{tmax:>10.4f}{tmax - filas[-1][1]:>+15.4f}{tt / 3600:>13.3f}{tf:>10.4f}{tg:>11.4f}{seg:>11.1f}")
        L.append("")
    L.append("  Las diferencias son indicadores de convergencia numerica, no cotas de error.")
    return L


def grafico_comparativo(ms, series):
    fig, (ax_q, ax_t) = plt.subplots(2, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    tt = np.linspace(0.0, max(DURACION_FUENTE_S, TIEMPO_FINAL_S) * 1.1, 600)
    ax_q.plot(tt / 3600, [fuente_q(float(x)) for x in tt], color="tab:red")
    ax_q.axvline(DURACION_FUENTE_S / 3600, color="gray", ls="--", label=f"tau = {DURACION_FUENTE_S / 3600:g} h")
    ax_q.axvline(TIEMPO_FINAL_S / 3600, color="k", ls=":", label=f"horizonte = {TIEMPO_FINAL_S / 3600:g} h")
    ax_q.set_ylabel("Q(t) [W/m³]"); ax_q.set_title(f"Fuente comun Q(t), Qmax = {QMAX_W_M3:.4f} W/m³"); ax_q.legend(); ax_q.grid(alpha=0.3)
    cmap = plt.get_cmap("tab10")
    for i, m in enumerate(ms):
        t, tc, tg = series[m["slug"]]
        ax_t.plot(t / 3600, tc, color=cmap(i), lw=2, label=m["titulo"])
        ax_t.plot(t / 3600, tg, color=cmap(i), lw=1, ls=":")
    ax_t.axhline(TEMPERATURA_REFERENCIA_C, color="tab:red", ls="--", label=f"referencia {TEMPERATURA_REFERENCIA_C:g} °C")
    ax_t.set_xlabel("Tiempo [h]"); ax_t.set_ylabel("Temperatura [°C]"); ax_t.grid(alpha=0.3); ax_t.legend(fontsize=8)
    ax_t.set_title("Temperatura central (trazo continuo) y maximo del dominio (punteado)")
    ruta = CARPETA_RESULTADOS / "comparacion" / "grafico_comparativo.png"
    fig.savefig(ruta, dpi=110)
    return fig, ruta


def crear_animacion(m, series, cortes, dt, escala, carpeta):
    """Tres cortes (xy, xz, yz) por el centro + curva central sincronizada; escala comun a todos los escenarios."""
    vmin, vmax = escala
    t_h = series[0] / 3600.0
    tc, tg = series[1], series[2]
    xy, xz, yz = cortes["xy"], cortes["xz"], cortes["yz"]
    tiempos_f, pasos_f = cortes["tiempos"], cortes["pasos"]
    fig = plt.figure(figsize=(13.5, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.9], width_ratios=[1.0, 1.15])
    ax_xy = fig.add_subplot(gs[0:2, 0]); ax_xz = fig.add_subplot(gs[0, 1]); ax_yz = fig.add_subplot(gs[1, 1]); ax_t = fig.add_subplot(gs[2, :])
    ims = []
    for ax, datos, ext, xl, yl, desc in ((ax_xy, xy, (0, LADO_X_M, 0, LADO_Y_M), "x [m]", "y [m]", cortes["desc"][0]),
                                         (ax_xz, xz, (0, LADO_X_M, 0, LADO_Z_M), "x [m]", "z [m]", cortes["desc"][1]),
                                         (ax_yz, yz, (0, LADO_Y_M, 0, LADO_Z_M), "y [m]", "z [m]", cortes["desc"][2])):
        im = ax.imshow(datos[0].T, origin="lower", extent=ext, vmin=vmin, vmax=vmax, cmap="inferno", aspect="equal")
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title("corte " + desc, fontsize=8)
        ims.append(im)
    fig.colorbar(ims[0], ax=[ax_xy, ax_xz, ax_yz], label="T [°C]", fraction=0.03, pad=0.02)
    ax_xy.plot(PUNTO_CENTRO_M[0], PUNTO_CENTRO_M[1], "o", ms=8, mfc="tab:cyan", mec="white")
    ax_xz.plot(PUNTO_CENTRO_M[0], PUNTO_CENTRO_M[2], "o", ms=6, mfc="tab:cyan", mec="white")
    ax_yz.plot(PUNTO_CENTRO_M[1], PUNTO_CENTRO_M[2], "o", ms=6, mfc="tab:cyan", mec="white")
    ax_t.set_xlim(0.0, t_h[-1]); ax_t.set_ylim(vmin - 1, vmax + 2)
    ax_t.axhline(TEMPERATURA_REFERENCIA_C, color="tab:red", ls="--", lw=1, label=f"referencia {TEMPERATURA_REFERENCIA_C:g} °C")
    l_c, = ax_t.plot([], [], color="tab:cyan", lw=2, label="centro")
    l_g, = ax_t.plot([], [], color="tab:orange", lw=1, ls=":", label="maximo del dominio")
    marca, = ax_t.plot([], [], "o", color="tab:cyan", ms=8, mec="white")
    ax_t.set_xlabel("Tiempo [h]"); ax_t.set_ylabel("T [°C]"); ax_t.grid(alpha=0.3); ax_t.legend(loc="upper left", fontsize=8)
    supt = fig.suptitle("")

    def actualizar(i):
        k = pasos_f[i]
        for im, datos in zip(ims, (xy, xz, yz)):
            im.set_data(datos[i].T)
        supt.set_text(f"{m['titulo']}: t = {tiempos_f[i] / 3600:.2f} h, Q = {fuente_q(float(tiempos_f[i])):.1f} W/m³, T centro = {tc[k]:.2f} °C")
        l_c.set_data(t_h[:k + 1], tc[:k + 1]); l_g.set_data(t_h[:k + 1], tg[:k + 1]); marca.set_data([t_h[k]], [tc[k]])
        return (*ims, supt, l_c, l_g, marca)

    actualizar(0)
    anim = FuncAnimation(fig, actualizar, frames=len(tiempos_f), interval=1000.0 / FPS_ANIMACION, blit=False, repeat=True)
    _ANIMACIONES.append(anim)
    if GUARDAR_ANIMACION:
        ruta = carpeta / f"evolucion_{m['slug']}.gif"
        anim.save(ruta, writer=PillowWriter(fps=FPS_ANIMACION), dpi=70)
        print(f"  GIF guardado: {ruta}")
    if not MOSTRAR_ANIMACION:
        plt.close(fig)
    return anim


def main():
    if NUM_FOTOGRAMAS < 2:
        raise ValueError("NUM_FOTOGRAMAS debe ser >= 2")
    CARPETA_RESULTADOS.mkdir(parents=True, exist_ok=True)
    escenarios = definir_escenarios()
    dt, pasos_frame, limites = paso_comun(escenarios)
    dx = (LADO_X_M / CELDAS_X, LADO_Y_M / CELDAS_Y, LADO_Z_M / CELDAS_Z)
    print(f"Caja {LADO_X_M:g}x{LADO_Y_M:g}x{LADO_Z_M:g} m, malla {CELDAS_X}x{CELDAS_Y}x{CELDAS_Z}; dt comun = {dt:.3f} s "
          f"(limites por escenario: " + ", ".join(f"{k} {v:.1f} s" for k, v in limites.items()) + ")")

    ms, series, cortes_por = [], {}, {}
    animados = [e["slug"] for e in escenarios] if ESCENARIOS_ANIMADOS is None else list(ESCENARIOS_ANIMADOS)
    textos = {}
    for esc in escenarios:
        sistema, r, seg = simular(esc, dt, capturar=True)
        m = metricas(esc, sistema, r, dt, seg)
        ms.append(m)
        series[esc["slug"]] = (r.historial_tiempos_s, r.historial_puntos_c[:, 0], r.historial_max_c)
        xy, xz, yz, descs = [], [], [], None
        for campo in r.estados_capturados_c:                 # solo los fotogramas; luego se descartan los campos 3D
            (a, b, c), descs = tres_cortes(campo, dx)
            xy.append(a); xz.append(b); yz.append(c)
        cortes_por[esc["slug"]] = dict(xy=np.array(xy), xz=np.array(xz), yz=np.array(yz), desc=descs,
                                       tiempos=np.array(r.instantes_capturados_s),
                                       pasos=[int(round(x / dt)) for x in r.instantes_capturados_s])
        carpeta = CARPETA_RESULTADOS / esc["slug"]
        carpeta.mkdir(parents=True, exist_ok=True)
        textos[esc["slug"]] = texto_escenario(m, dt, pasos_frame, descs)
        (carpeta / "resumen.txt").write_text(textos[esc["slug"]] + "\n", encoding="utf-8")
        guardar_csv_escenario(carpeta, r)
        print(f"  [{esc['slug']}] Tc max {m['Tc_max']:.3f} °C en {m['t_Tc_max_s'] / 3600:.2f} h; final {m['Tc_fin']:.3f} °C; {seg:.1f} s")
        del r

    # Verificaciones y sensibilidad
    verif = verificaciones_rapidas()
    sens = sensibilidad(escenarios, dt) if EJECUTAR_SENSIBILIDAD else None
    neu = next(m for m in ms if m["slug"] == "neumann_cero")

    L = ["EJERCICIO FUENTE VARIABLE 3D - COMPARACION DE ESCENARIOS", "",
         "PARAMETROS",
         f"  Caja {LADO_X_M:g} x {LADO_Y_M:g} x {LADO_Z_M:g} m; malla {CELDAS_X}x{CELDAS_Y}x{CELDAS_Z} (dx, dy, dz = {dx[0]:g}, {dx[1]:g}, {dx[2]:g} m); centro {PUNTO_CENTRO_M} m",
         f"  T inicial {TEMPERATURA_INICIAL_C:g} °C en todos los escenarios; referencia {TEMPERATURA_REFERENCIA_C:g} °C (superarla se INFORMA, Qmax no se reajusta)",
         f"  k = {CONDUCTIVIDAD_W_M_K:g} W/(m K), rho = {DENSIDAD_KG_M3:g} kg/m³, cp = {CALOR_ESPECIFICO_J_KG_K:g} J/(kg K) (los del ejercicio 2D)",
         f"  Qmax = {QMAX_W_M3!r} W/m³ (recuperado del ejercicio Dirichlet 2D; no recalculado); tau = {DURACION_FUENTE_S:g} s; horizonte = {TIEMPO_FINAL_S / 3600:g} h",
         f"  dt efectivo comun = {dt:.3f} s ({pasos_frame} pasos por fotograma, {NUM_FOTOGRAMAS} fotogramas, FPS solo reproduccion = {FPS_ANIMACION}); limites: " + ", ".join(f"{k} {v:.1f} s" for k, v in limites.items()),
         f"  Robin: beta = h/k [1/m] = {BETAS_ROBIN_1_M} (EXPLORATORIOS, no datos del profesor); h = k·beta; -k dT/dn_ext = h (T_s - T_amb), ambiente {TEMPERATURA_AMBIENTE_C:g} °C",
         "", "TABLA COMPARATIVA"] + tabla_texto(ms) + [
         "", "BALANCE DE ENERGIA POR ESCENARIO (J; frontera positiva = entra al cuerpo, negativa = sale)",
         f"  {'escenario':<18}{'cambio':>14}{'fuente':>14}{'fronteras':>14}{'residuo':>12}"]
    for m in ms:
        L.append(f"  {m['slug']:<18}{m['e_fin'] - m['e_ini']:>+14.4e}{m['e_fuente']:>+14.4e}{m['e_front']:>+14.4e}{m['residuo']:>+12.2e}")
    err_neu = float(np.max(np.abs(series["neumann_cero"][1] - temperatura_adiabatica_exacta(series["neumann_cero"][0]))))
    L += ["", "VERIFICACIONES",
          f"  Neumann cero (malla completa): T final {neu['Tc_fin']:.4f} °C (esperado ~71.30 °C); dispersion {neu['dispersion_final']:.2e} °C; "
          f"error maximo vs solucion adiabatica {err_neu:.3e} °C (tolerancia {tolerancia_integracion_c(dt):.3e} °C) -> {'CUMPLE' if err_neu <= tolerancia_integracion_c(dt) else 'NO CUMPLE'}"]
    for nombre, ok, det in verif:
        L.append(f"  {nombre}: {'CUMPLE' if ok else 'NO CUMPLE'} ({det}) [malla {'x'.join(map(str, CELDAS_VERIFICACION))}]")
    L.append(f"  Simetria final (peor escenario): {max(m['simetria'] for m in ms):.2e} °C")
    if sens:
        L += ["", "SENSIBILIDAD (metricas centrales; costo = tiempo de calculo)"] + sens
    texto = "\n".join(L)
    (CARPETA_RESULTADOS / "comparacion").mkdir(parents=True, exist_ok=True)
    (CARPETA_RESULTADOS / "comparacion" / "resumen.txt").write_text(texto + "\n", encoding="utf-8")
    print(texto)
    print(f"Resumen guardado: {CARPETA_RESULTADOS / 'comparacion' / 'resumen.txt'}")
    print(f"Tabla comparativa: {guardar_tabla_comparativa(ms)}")

    fig, ruta = grafico_comparativo(ms, series)
    print(f"Grafico comparativo: {ruta}")
    if not MOSTRAR_GRAFICOS:
        plt.close(fig)

    if MOSTRAR_ANIMACION or GUARDAR_ANIMACION:
        # Escala de color COMUN a todos los escenarios (y que incluye la referencia de 70 °C).
        todos = [c[k].min() for c in cortes_por.values() for k in ("xy", "xz", "yz")] + [c[k].max() for c in cortes_por.values() for k in ("xy", "xz", "yz")]
        escala = (float(min(todos + [TEMPERATURA_INICIAL_C])), float(max(todos + [TEMPERATURA_REFERENCIA_C])))
        print(f"Escala de color comun: {escala[0]:.2f} a {escala[1]:.2f} °C")
        for m in ms:
            if m["slug"] in animados:
                crear_animacion(m, series[m["slug"]], cortes_por[m["slug"]], dt, escala, CARPETA_RESULTADOS / m["slug"])
    if MOSTRAR_GRAFICOS or MOSTRAR_ANIMACION:
        plt.show()      # una sola llamada para todas las ventanas
    if not (all(ok for _, ok, _ in verif) and err_neu <= tolerancia_integracion_c(dt)):
        print("VERIFICACION FALLIDA: revisa el resumen.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
