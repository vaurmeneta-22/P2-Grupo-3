"""Análisis de resultados 3D para el entregable E09.

Las funciones de este módulo no simulan ni modifican el solver: solo extraen
series temporales, seleccionan puntos representativos, calculan cortes y
resumen métricas a partir de un :class:`ResultadoEtapas` ya calculado.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mcplan.geometry.voxelizacion import GrillaVoxel, caras_superficie_exterior
from .vertidas import ResultadoEtapas


@dataclass(frozen=True)
class PuntoAnalisis:
    """Punto de la grilla elegido para una curva de temperatura representativa."""

    etiqueta: str
    categoria: str
    ijk: tuple[int, int, int]
    xyz_m: tuple[float, float, float]
    vertida_id: int


@dataclass(frozen=True)
class CorteTermico:
    """Corte 2D de un campo térmico en un plano de la grilla."""

    eje: str
    indice: int
    coordenada_m: float
    tiempo_s: float
    tiempo_h: float
    temperaturas_c: np.ndarray
    mascara_activa: np.ndarray
    extension_m: tuple[float, float]


@dataclass(frozen=True)
class MetricasTmax:
    """Resumen global de la serie de temperatura máxima."""

    t_pico_s: float
    t_pico_h: float
    tmax_pico_c: float
    tmax_final_c: float
    tmedia_final_c: float
    n_muestras: int


def _validar_ijk(ijk: tuple[int, int, int], forma: tuple[int, int, int]) -> None:
    if len(ijk) != 3:
        raise ValueError("ijk debe contener tres índices")
    for indice, limite in zip(ijk, forma):
        if not isinstance(indice, (int, np.integer)) or not 0 <= int(indice) < limite:
            raise ValueError(f"ijk={ijk} fuera de la grilla de forma {forma}")


def serie_punto(resultado: ResultadoEtapas,
                ijk: tuple[int, int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Temperatura [°C] en un voxel para cada instante guardado.

    Devuelve ``(horas, temperaturas)``. La temperatura es ``NaN`` mientras el
    voxel no pertenece a una vertida activa, de modo que la curva no inventa
    valores antes de la colocación.
    """
    if not resultado.temperaturas_c:
        raise ValueError("El resultado no contiene estados térmicos")
    forma = resultado.temperaturas_c[0].shape
    _validar_ijk(ijk, forma)
    horas = np.asarray(resultado.tiempos_s, dtype=float) / 3600.0
    temperaturas = np.full(horas.shape, np.nan, dtype=float)
    for indice, (campo, activa) in enumerate(
            zip(resultado.temperaturas_c, resultado.mascaras_activas)):
        if bool(activa[ijk]):
            temperaturas[indice] = float(campo[ijk])
    return horas, temperaturas


def _voxel_mas_caliente(campo: np.ndarray, activa: np.ndarray) -> tuple[int, int, int]:
    if not np.any(activa):
        raise ValueError("No hay voxeles activos para buscar el máximo")
    valores = np.where(activa, campo, -np.inf)
    indices = np.unravel_index(int(np.argmax(valores)), campo.shape)
    return tuple(int(indice) for indice in indices)


def _picos_por_voxel(resultado: ResultadoEtapas,
                     mascara: np.ndarray) -> np.ndarray:
    picos = np.full(mascara.shape, -np.inf, dtype=float)
    for campo, activa in zip(resultado.temperaturas_c, resultado.mascaras_activas):
        seleccion = mascara & activa
        if np.any(seleccion):
            picos[seleccion] = np.maximum(picos[seleccion], campo[seleccion])
    return picos


def seleccionar_puntos_representativos(
        grilla: GrillaVoxel,
        resultado: ResultadoEtapas,
) -> tuple[PuntoAnalisis, ...]:
    """Elige tres puntos de forma determinista y trazable.

    * P1: voxel con la mayor temperatura de toda la simulación.
    * P2: voxel de superficie expuesta con la mayor temperatura de su historia.
    * P3: voxel interior de la última vertida activada, cercano a su centroide.

    La selección no depende de posiciones codificadas a mano y puede
    reproducirse con los mismos CSV de entrada.
    """
    if not resultado.temperaturas_c:
        raise ValueError("El resultado no contiene estados térmicos")
    if len(resultado.temperaturas_c) != len(resultado.mascaras_activas):
        raise ValueError("Estados y máscaras activas deben tener la misma longitud")

    indice_pico = int(np.argmax(resultado.temperaturas_maximas_c))
    ijk_p1 = _voxel_mas_caliente(
        resultado.temperaturas_c[indice_pico],
        resultado.mascaras_activas[indice_pico],
    )

    activa_final = resultado.mascaras_activas[-1]
    caras = caras_superficie_exterior(grilla.vertida_id, activa_final)
    superficie = np.logical_or.reduce(tuple(caras.values()))
    picos_superficie = _picos_por_voxel(resultado, superficie)
    candidatos_superficie = np.isfinite(picos_superficie)
    if np.any(candidatos_superficie):
        # Evita repetir P1 cuando exista otra cara expuesta con historia térmica.
        sin_p1 = candidatos_superficie.copy()
        sin_p1[ijk_p1] = False
        if np.any(sin_p1):
            candidatos_superficie = sin_p1
        ijk_p2 = tuple(int(i) for i in np.unravel_index(
            int(np.argmax(np.where(candidatos_superficie, picos_superficie, -np.inf))),
            picos_superficie.shape,
        ))
    else:
        ijk_p2 = ijk_p1

    ids_activos = [int(i) for i in np.unique(grilla.vertida_id[activa_final]) if i > 0]
    if not ids_activos:
        raise ValueError("No hay vertidas activas en el estado final")
    id_ultima = max(ids_activos)
    pertenece_ultima = grilla.vertida_id == id_ultima
    interior = pertenece_ultima & activa_final & ~superficie
    seleccion = interior if np.any(interior) else (pertenece_ultima & activa_final)
    posiciones = np.argwhere(seleccion)
    centroide = posiciones.mean(axis=0)
    distancias = np.sum((posiciones - centroide) ** 2, axis=1)
    ijk_p3 = tuple(int(i) for i in posiciones[int(np.argmin(distancias))])

    puntos = []
    for etiqueta, categoria, ijk in (
            ("P1", "pico_global", ijk_p1),
            ("P2", "superficie_expuesta", ijk_p2),
            ("P3", "nucleo_ultima_vertida", ijk_p3),
    ):
        xyz = tuple(
            float(grilla.origen_m[d] + (ijk[d] + 0.5) * grilla.paso_m[d])
            for d in range(3)
        )
        puntos.append(PuntoAnalisis(
            etiqueta=etiqueta,
            categoria=categoria,
            ijk=ijk,
            xyz_m=xyz,
            vertida_id=int(grilla.vertida_id[ijk]),
        ))
    return tuple(puntos)


def corte_termico(grilla: GrillaVoxel,
                  resultado: ResultadoEtapas,
                  tiempo_s: float,
                  eje: str,
                  indice: int | None = None) -> CorteTermico:
    """Extrae un corte 2D en el instante guardado más cercano a ``tiempo_s``.

    ``eje`` es ``"x"``, ``"y"`` o ``"z"``. Si ``indice`` es ``None`` se usa un
    plano central que atraviesa el sólido activo. Las celdas inactivas quedan
    como ``NaN`` para no confundirlas con temperaturas frías.
    """
    eje = eje.lower()
    if eje not in ("x", "y", "z"):
        raise ValueError("eje debe ser 'x', 'y' o 'z'")
    if not np.isfinite(tiempo_s):
        raise ValueError("tiempo_s debe ser finito")
    if not resultado.tiempos_s.size:
        raise ValueError("El resultado no contiene tiempos")

    eje_indice = "xyz".index(eje)
    indice_tiempo = int(np.argmin(np.abs(resultado.tiempos_s - float(tiempo_s))))
    activa = resultado.mascaras_activas[indice_tiempo]
    campo = resultado.temperaturas_c[indice_tiempo]

    if indice is None:
        otros_ejes = tuple(d for d in range(3) if d != eje_indice)
        con_solido = np.any(activa, axis=otros_ejes)
        disponibles = np.flatnonzero(con_solido)
        if not disponibles.size:
            raise ValueError("No hay sólido activo en el instante solicitado")
        indice = int(disponibles[len(disponibles) // 2])
    indice = int(indice)
    if not 0 <= indice < grilla.forma[eje_indice]:
        raise ValueError(f"indice={indice} fuera del eje {eje}")

    rebanador = [slice(None)] * 3
    rebanador[eje_indice] = indice
    rebanador = tuple(rebanador)
    temperaturas = np.where(activa[rebanador], campo[rebanador], np.nan)

    otros = [d for d in range(3) if d != eje_indice]
    extension = tuple(
        float(grilla.origen_m[d] + grilla.forma[d] * grilla.paso_m[d] - grilla.origen_m[d])
        for d in otros
    )
    return CorteTermico(
        eje=eje,
        indice=indice,
        coordenada_m=float(grilla.origen_m[eje_indice]
                           + (indice + 0.5) * grilla.paso_m[eje_indice]),
        tiempo_s=float(resultado.tiempos_s[indice_tiempo]),
        tiempo_h=float(resultado.tiempos_s[indice_tiempo] / 3600.0),
        temperaturas_c=temperaturas,
        mascara_activa=np.asarray(activa[rebanador], dtype=bool).copy(),
        extension_m=extension,
    )


def cortes_por_punto(grilla: GrillaVoxel,
                     resultado: ResultadoEtapas,
                     ijk: tuple[int, int, int],
                     tiempo_s: float) -> tuple[CorteTermico, CorteTermico, CorteTermico]:
    """Tres cortes ortogonales que pasan por ``ijk`` en un mismo instante."""
    _validar_ijk(ijk, grilla.forma)
    return tuple(
        corte_termico(grilla, resultado, tiempo_s, eje, indice)
        for eje, indice in zip("xyz", ijk)
    )


def metricas_tmax(resultado: ResultadoEtapas) -> MetricasTmax:
    """Pico y valores finales de la serie ``Tmax(t)`` del solver."""
    maximas = np.asarray(resultado.temperaturas_maximas_c, dtype=float)
    if maximas.size == 0:
        raise ValueError("El resultado no contiene temperaturas máximas")
    if not np.any(np.isfinite(maximas)):
        raise ValueError("La serie de temperaturas máximas no contiene valores finitos")
    indice = int(np.nanargmax(maximas))
    return MetricasTmax(
        t_pico_s=float(resultado.tiempos_s[indice]),
        t_pico_h=float(resultado.tiempos_s[indice] / 3600.0),
        tmax_pico_c=float(maximas[indice]),
        tmax_final_c=float(maximas[-1]),
        tmedia_final_c=float(resultado.temperaturas_promedio_c[-1]),
        n_muestras=int(maximas.size),
    )


def eventos_activacion(resultado: ResultadoEtapas) -> tuple[tuple[float, tuple[int, ...]], ...]:
    """Eventos ``(hora, IDs)`` en los que se incorporó al menos una vertida."""
    return tuple(
        (float(tiempo / 3600.0), tuple(int(i) for i in ids))
        for tiempo, ids in zip(resultado.tiempos_s, resultado.ids_activados)
        if ids
    )
