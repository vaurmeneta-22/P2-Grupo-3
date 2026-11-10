"""Simulación 3D por etapas usando la grilla IFC y el solver multidimensional."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from .diffusion_nd import CondicionBorde, DifusionCartesiana, aislado
from .hidratacion import ModeloHidratacionExponencial, ModeloPulsoCosenoidal


@dataclass(frozen=True)
class ResultadoEtapas:
    tiempos_s: np.ndarray
    temperaturas_c: tuple[np.ndarray, ...]
    mascaras_activas: tuple[np.ndarray, ...]
    ids_activados: tuple[tuple[int, ...], ...]
    temperaturas_maximas_c: np.ndarray
    temperaturas_promedio_c: np.ndarray
    energia_fuente_j: float
    energia_activacion_j: float
    energia_frontera_j: float
    energia_final_j: float
    error_balance_j: float
    ultimo_paso_maximo_estable_s: float
    retiros: tuple = ()   # (t_s, nombres, dt_estable_antes_s, dt_estable_despues_s)


def simular_vertidas(*, grilla: GrillaVoxel, vertidas: list[VertidaIFC],
                     conductividad_w_m_k: float | np.ndarray,
                     densidad_kg_m3: float | np.ndarray,
                     calor_especifico_j_kg_k: float | np.ndarray,
                     hidratacion: ModeloHidratacionExponencial | ModeloPulsoCosenoidal,
                     tiempo_final_s: float, fraccion_paso_estable: float = 0.8,
                     paso_tiempo_max_s: float = 600.0,
                     tiempos_muestreo_s: tuple[float, ...] | list[float] | np.ndarray | None = None,
                     fronteras: dict[str, CondicionBorde] | None = None,
                     temperatura_inicial_c: float = 20.0,
                     progreso: Callable[[str], None] | None = None,
                     seguimiento=None, bordes=None) -> ResultadoEtapas:
    """Activa cada vertida en su hora IFC, la inicializa y simula conducción.

    La lista de vertidas constituye la entrada del plan. El solver recibe una
    máscara activa que cambia en los eventos, y no contiene horas de inicio
    codificadas. Las superficies hacia voxeles inactivos usan ``fronteras``;
    por defecto son adiabáticas. Las interfaces entre dos celdas activas
    intercambian calor con la ley de conductancia del solver existente.

    ``seguimiento`` (opcional: un observador o una lista de ellos, ver
    ``seguimiento.SeguimientoPuntos`` y ``extremos.SeguimientoExtremos``) solo observa
    el campo tras cada paso y en cada activación; no modifica la solución.

    ``bordes`` (opcional, ``bordes.ModeloBordes``) aporta suelo en la base y
    resistencias de protección (moldaje, aislación) por vertida y cara. Cada retiro es
    un instante de evento: se llega a él integrando, se actualizan las conductancias
    conservando el campo de temperatura y el siguiente intervalo recalcula el paso
    estable. Con ``bordes=None`` el comportamiento es el anterior.
    """
    if tiempo_final_s <= 0 or not np.isfinite(tiempo_final_s):
        raise ValueError("tiempo_final_s debe ser positivo y finito")
    if not 0 < fraccion_paso_estable <= 1:
        raise ValueError("fraccion_paso_estable debe pertenecer a (0, 1]")
    if paso_tiempo_max_s <= 0 or not np.isfinite(paso_tiempo_max_s):
        raise ValueError("paso_tiempo_max_s debe ser positivo y finito")
    if not vertidas:
        raise ValueError("Se requiere al menos una vertida")
    if tiempos_muestreo_s is None:
        muestras_s = np.empty(0, dtype=float)
    else:
        muestras_s = np.asarray(tiempos_muestreo_s, dtype=float)
        if muestras_s.ndim != 1 or not np.all(np.isfinite(muestras_s)):
            raise ValueError("tiempos_muestreo_s debe ser una lista unidimensional finita")
        if np.any(muestras_s < 0.0) or np.any(muestras_s > tiempo_final_s):
            raise ValueError("Los tiempos de muestreo deben estar entre 0 y tiempo_final_s")

    ids = grilla.vertida_id
    forma = grilla.forma
    inicio_s = np.full(forma, np.inf)
    colocacion_c = np.full(forma, float(temperatura_inicial_c))
    por_id = {v.vertida_id: v for v in vertidas}
    if set(grilla.volumen_ifc_por_vertida_m3) - set(por_id):
        raise ValueError("Faltan metadatos de una o más vertidas presentes en la grilla")
    for identificador in np.unique(ids):
        if identificador == 0:
            continue
        if int(identificador) not in por_id:
            raise ValueError(f"No hay metadatos para VertidaID={identificador}")
        vertida = por_id[int(identificador)]
        seleccion = ids == identificador
        inicio_s[seleccion] = vertida.inicio_h * 3600.0
        colocacion_c[seleccion] = vertida.temperatura_colocacion_c

    rho = np.broadcast_to(np.asarray(densidad_kg_m3, dtype=float), forma)
    cp = np.broadcast_to(np.asarray(calor_especifico_j_kg_k, dtype=float), forma)
    if np.any(rho <= 0) or np.any(cp <= 0):
        raise ValueError("densidad y calor específico deben ser positivos")
    energia_volumetrica_j_k = rho * cp * grilla.volumen_voxel_m3

    def qdot(tiempo_s):
        edad_s = tiempo_s - inicio_s
        return hidratacion.generacion_w_m3(edad_s, rho, cp)

    mascara = (ids > 0) & (inicio_s <= 0.0)
    temperatura = np.zeros(forma, dtype=float)
    temperatura[mascara] = colocacion_c[mascara]
    if bordes is not None:
        condiciones = bordes.fronteras()
    else:
        condiciones = fronteras or {cara: aislado() for cara in ("x-", "x+", "y-", "y+", "z-", "z+")}
    solver = DifusionCartesiana(
        forma=forma,
        dx_m=grilla.paso_m,
        conductividad_w_m_k=conductividad_w_m_k,
        densidad_kg_m3=rho,
        calor_especifico_j_kg_k=cp,
        fuente_w_m3=qdot,
        fronteras=condiciones,
        mascara_activa=mascara,
        resistencias_borde=None if bordes is None else bordes.resistencias_en(0.0),
        fronteras_celda=None if bordes is None else bordes.fronteras_celda(),
    )

    tiempos_guardados = [0.0]
    temperaturas_guardadas = [temperatura.copy()]
    mascaras_guardadas = [mascara.copy()]
    ids_evento = [tuple(sorted(int(i) for i in np.unique(ids[mascara]) if i > 0))]
    maxima, promedios = [], []
    energia_activacion = float(np.sum(energia_volumetrica_j_k[mascara]
                                      * temperatura[mascara]))
    energia_fuente = 0.0
    energia_frontera = 0.0
    ultimo_dt_estable = np.inf

    def guardar_estado(tiempo, nuevos=()):
        activos = solver.mascara_activa
        if np.any(activos):
            maxima.append(float(np.max(temperatura[activos])))
            promedios.append(float(np.mean(temperatura[activos])))
        else:
            maxima.append(float("nan"))
            promedios.append(float("nan"))
        if tiempos_guardados[-1] != tiempo:
            tiempos_guardados.append(float(tiempo))
            temperaturas_guardadas.append(temperatura.copy())
            mascaras_guardadas.append(activos.copy())
            ids_evento.append(tuple(sorted(int(i) for i in nuevos)))
        else:
            temperaturas_guardadas[-1] = temperatura.copy()
            mascaras_guardadas[-1] = activos.copy()
            ids_evento[-1] = tuple(sorted(int(i) for i in nuevos))

    guardar_estado(0.0, ids_evento[0])
    seguimientos = ([] if seguimiento is None
                    else list(seguimiento) if isinstance(seguimiento, (list, tuple))
                    else [seguimiento])

    def registrar_seguimientos(tiempo_s, campo):
        for observado in seguimientos:
            observado.registrar(tiempo_s, campo, solver)

    registrar_seguimientos(0.0, temperatura)
    observador = registrar_seguimientos if seguimientos else None
    eventos_h = sorted({float(v.inicio_h) for v in vertidas if 0.0 < v.inicio_h * 3600.0 < tiempo_final_s})
    eventos_s = [h * 3600.0 for h in eventos_h] + [float(tiempo_final_s)]
    retiros_s = [] if bordes is None else bordes.instantes_retiro_s(tiempo_final_s)
    eventos_s += retiros_s
    retiros_registrados = []
    # El estado inicial ya se guarda antes del ciclo; evitar duplicar t=0.
    instantes_guardado_s = sorted({tiempo for tiempo in eventos_s + muestras_s.tolist()
                                   if tiempo > 0.0})
    tiempo_actual = 0.0
    identificadores_activados = set(int(i) for i in np.unique(ids[mascara]) if i > 0)

    for tiempo_evento in instantes_guardado_s:
        intervalo = tiempo_evento - tiempo_actual
        if intervalo > 0 and np.any(solver.mascara_activa):
            ultimo_dt_estable = solver.paso_maximo_estable_s()
            dt_estable = (intervalo if not np.isfinite(ultimo_dt_estable)
                          else fraccion_paso_estable * ultimo_dt_estable)
            dt_max = min(float(paso_tiempo_max_s), dt_estable)
            resultado = solver.simular(
                temperatura_inicial_c=temperatura,
                tiempo_final_s=intervalo,
                paso_tiempo_max_s=dt_max,
                tiempo_inicial_s=tiempo_actual,
                observador=observador,
            )
            temperatura = resultado.temperatura_c
            energia_fuente += resultado.energia_fuente_j
            energia_frontera += resultado.energia_frontera_total_j
        tiempo_actual = tiempo_evento

        nuevos = []
        for vertida in vertidas:
            if vertida.inicio_h * 3600.0 <= tiempo_actual + 1e-9 and vertida.vertida_id not in identificadores_activados:
                celdas_nuevas = ids == vertida.vertida_id
                temperatura[celdas_nuevas] = vertida.temperatura_colocacion_c
                energia_activacion += float(np.sum(energia_volumetrica_j_k[celdas_nuevas]
                                                   * temperatura[celdas_nuevas]))
                identificadores_activados.add(vertida.vertida_id)
                nuevos.append(vertida.vertida_id)
        mascara = (ids > 0) & (inicio_s <= tiempo_actual + 1e-9)
        solver.actualizar_mascara_activa(mascara)
        if bordes is not None:
            retiradas = bordes.retiradas_en(tiempo_actual)
            dt_antes = solver.paso_maximo_estable_s()
            bordes.aplicar(solver, tiempo_actual)      # T no cambia; cambian las conductancias
            if retiradas:
                dt_despues = solver.paso_maximo_estable_s()
                retiros_registrados.append((float(tiempo_actual),
                                            tuple(f"V{p.vertida_id}:{p.nombre}" for p in retiradas),
                                            float(dt_antes), float(dt_despues)))
                if progreso:
                    progreso(f"t={tiempo_actual / 3600:.2f} h; retiro de "
                             f"{[f'V{p.vertida_id}:{p.nombre}' for p in retiradas]}; paso estable "
                             f"{dt_antes:.1f} s -> {dt_despues:.1f} s")
        guardar_estado(tiempo_actual, nuevos)
        # Reemplaza la muestra de este instante: ahora incluye vertidas recién activadas.
        registrar_seguimientos(tiempo_actual, temperatura)
        if progreso:
            progreso(f"t={tiempo_actual / 3600:.2f} h; activadas={nuevos}; "
                     f"voxeles activos={int(mascara.sum())}")

    energia_final = solver.energia_j(temperatura)
    balance = energia_final - energia_activacion - energia_fuente - energia_frontera
    return ResultadoEtapas(
        tiempos_s=np.asarray(tiempos_guardados),
        temperaturas_c=tuple(temperaturas_guardadas),
        mascaras_activas=tuple(mascaras_guardadas),
        ids_activados=tuple(ids_evento),
        temperaturas_maximas_c=np.asarray(maxima),
        temperaturas_promedio_c=np.asarray(promedios),
        energia_fuente_j=energia_fuente,
        energia_activacion_j=energia_activacion,
        energia_frontera_j=energia_frontera,
        energia_final_j=energia_final,
        error_balance_j=balance,
        ultimo_paso_maximo_estable_s=float(ultimo_dt_estable),
        retiros=tuple(retiros_registrados),
    )
