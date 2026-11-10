"""Seguimiento de temperaturas puntuales y de la diferencia interior-superficie.

Para cada vertida se fijan DOS puntos durante toda la corrida:

* un punto interior (centro de un voxel que pertenece a la vertida), y
* un voxel próximo a la superficie: el voxel de borde de esa vertida junto a una
  cara ``cara`` (``x-``, ``x+``, ``y-``, ``y+``, ``z-`` o ``z+``). Su CENTRO es
  distinto del punto de la CARA exterior, donde se reconstruye ``T_cara_estimada``
  (a d/2 del centro, sobre la normal de la cara).

Los puntos se dan a mano (``resolver_puntos``, validación estricta) o se proponen
para cada grilla (``proponer_puntos``, modo automático). En ambos casos quedan
fijos durante una corrida; en modo automático sus posiciones físicas cambian al
cambiar dx, así que las diferencias entre resoluciones NO comparan el mismo punto.

El seguimiento NO integra nada: solo observa el campo que ya calcula el solver
en cada paso de integración. No se guardan campos completos adicionales ni se
ejecuta una simulación por punto.

TEMPERATURA SUPERFICIAL (dos magnitudes distintas)
--------------------------------------------------
* ``T_celda_proxima_superficie``: temperatura del CENTRO del voxel de borde.
  Está a d/2 de la cara exterior (d = ancho del voxel en la dirección de la cara).
* ``T_cara_estimada``: temperatura reconstruida EN la cara exterior con las mismas
  resistencias de media celda que usa el solver, ``Rc = d/(2 k)``. Con ``q_in`` el
  flujo entrante [W/m²] por esa cara:

      Robin:     q_in = (T_amb - T_celda)/(Rc + 1/h)   ->  T_cara = T_celda + q_in Rc
      Dirichlet: T_cara = T prescrita
      Neumann:   q_in dado                              ->  T_cara = T_celda + q_in Rc

  Es una reconstrucción 1D normal a la cara, consistente con el esquema, no una
  medición ni una solución exacta de la temperatura superficial.

Si la cara deja de ser exterior porque se activa otra vertida al lado, el estado
pasa a ``interfaz con vertida N`` y ``T_cara_estimada`` queda en n/a: esa cara ya
no está expuesta. Antes de activarse la vertida todo es n/a (nunca 0 °C).

La diferencia es ``dT = T_interior - T_superficie`` entre dos puntos [°C o K],
no un gradiente en K/m. Es una métrica térmica, no un análisis mecánico.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np

from .diffusion_nd import CondicionBorde

CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")
NA = "n/a"


@dataclass(frozen=True)
class PuntoVertida:
    """Par de puntos fijos (interior, superficie) de una vertida, ya validados."""

    vertida_id: int
    indice_interior: tuple[int, int, int]
    indice_voxel_proximo: tuple[int, int, int]
    cara: str
    interior_solicitado_m: tuple[float, float, float]
    voxel_proximo_solicitado_m: tuple[float, float, float]
    centro_interior_m: tuple[float, float, float]
    centro_voxel_proximo_m: tuple[float, float, float]
    centro_cara_m: tuple[float, float, float]
    vecino_id: int  # vertida al otro lado de la cara (0 = vacío o borde de la grilla)

    @property
    def distancia_interior_a_cara_m(self) -> float:
        """Distancia entre el centro del voxel interior y el centro de la cara medida."""
        return float(np.linalg.norm(np.subtract(self.centro_cara_m, self.centro_interior_m)))


def _indice_de_punto(grilla, punto_m, nombre):
    punto = np.asarray(punto_m, dtype=float)
    if punto.shape != (3,) or not np.all(np.isfinite(punto)):
        raise ValueError(f"{nombre}: se esperan tres coordenadas finitas [x, y, z] en metros")
    indice = np.floor((punto - np.asarray(grilla.origen_m)) / np.asarray(grilla.paso_m)).astype(int)
    if np.any(indice < 0) or np.any(indice >= np.asarray(grilla.forma)):
        raise ValueError(f"{nombre}: el punto {punto.tolist()} cae fuera de la grilla")
    return tuple(int(i) for i in indice)


def _centro_voxel(grilla, indice):
    return tuple(float(grilla.origen_m[d] + (indice[d] + 0.5) * grilla.paso_m[d])
                 for d in range(3))


def _vecino(grilla, indice, cara):
    eje, signo = "xyz".index(cara[0]), (-1 if cara[1] == "-" else 1)
    vecino = list(indice)
    vecino[eje] += signo
    if vecino[eje] < 0 or vecino[eje] >= grilla.forma[eje]:
        return 0
    return int(grilla.vertida_id[tuple(vecino)])


def resolver_puntos(grilla, puntos_por_vertida: dict) -> list[PuntoVertida]:
    """Convierte la configuración en puntos validados sobre la grilla.

    Valida que ambos puntos pertenezcan a su vertida y que la cara indicada no
    sea interna a ella. Un vecino de OTRA vertida no es error: la cara empieza
    expuesta y pasa a interfaz cuando esa vertida se activa (el seguimiento lo
    registra); ``vecino_id`` lo informa.
    """
    resultado = []
    for clave, cfg in sorted(puntos_por_vertida.items(), key=lambda par: int(par[0])):
        identificador = int(clave)
        for campo in ("interior_m", "voxel_proximo_m", "cara_superficie"):
            if campo not in cfg:
                raise ValueError(f"Vertida {identificador}: falta '{campo}' en seguimiento_puntos")
        cara = str(cfg["cara_superficie"])
        if cara not in CARAS:
            raise ValueError(f"Vertida {identificador}: cara_superficie debe ser una de {CARAS}")
        interior = _indice_de_punto(grilla, cfg["interior_m"], f"Vertida {identificador} interior_m")
        superficie = _indice_de_punto(grilla, cfg["voxel_proximo_m"], f"Vertida {identificador} voxel_proximo_m")
        for nombre, indice in (("interior_m", interior), ("voxel_proximo_m", superficie)):
            if int(grilla.vertida_id[indice]) != identificador:
                raise ValueError(
                    f"Vertida {identificador}: {nombre} cae en el voxel {indice}, que pertenece a "
                    f"VertidaID={int(grilla.vertida_id[indice])}; ajusta la coordenada o usa modo 'automatico'")
        vecino = _vecino(grilla, superficie, cara)
        if vecino == identificador:
            raise ValueError(f"Vertida {identificador}: la cara {cara} del voxel {superficie} es "
                             "interna a la vertida, no una superficie")
        centro_voxel_proximo = _centro_voxel(grilla, superficie)
        eje, signo = "xyz".index(cara[0]), (-1 if cara[1] == "-" else 1)
        centro_cara = list(centro_voxel_proximo)
        centro_cara[eje] += signo * 0.5 * grilla.paso_m[eje]
        resultado.append(PuntoVertida(
            vertida_id=identificador, indice_interior=interior, indice_voxel_proximo=superficie,
            cara=cara,
            interior_solicitado_m=tuple(float(v) for v in cfg["interior_m"]),
            voxel_proximo_solicitado_m=tuple(float(v) for v in cfg["voxel_proximo_m"]),
            centro_interior_m=_centro_voxel(grilla, interior),
            centro_voxel_proximo_m=centro_voxel_proximo,
            centro_cara_m=tuple(centro_cara), vecino_id=vecino))
    return resultado


CARAS_LATERALES = ("x+", "x-", "y+", "y-")


def proponer_puntos(grilla, cara: str | None = "x+") -> dict:
    """Propone puntos a partir de los voxeles reales de cada vertida (modo automático).

    Interior: voxel de la vertida más cercano al eje central (centro de su caja
    de voxeles en x e y) a media altura de la vertida. Voxel próximo a la
    superficie: avanzando desde el interior en la dirección de ``cara`` (misma
    altura y mismo eje), el último voxel de la vertida. Se exige que la cara quede
    libre (sin voxel de ninguna vertida al otro lado); si ``cara`` es None se
    prueba x+, x-, y+ e y- y se usa la primera que cumpla, y si se fuerza una
    ``cara`` que no cumple se lanza ValueError. Devuelve el diccionario
    ``puntos_por_vertida`` en coordenadas globales [m] (centros de voxel).

    Los puntos dependen de dx: al cambiar la resolución la propuesta cae en otros
    voxeles y, en general, en otras posiciones físicas.
    """
    if cara is not None and (cara not in CARAS or cara[0] == "z"):
        raise ValueError("la cara propuesta debe ser lateral (x-, x+, y- o y+) o None")
    candidatas = CARAS_LATERALES if cara is None else (cara,)
    propuesta = {}
    for identificador in sorted(int(i) for i in np.unique(grilla.vertida_id) if i > 0):
        indices = np.argwhere(grilla.vertida_id == identificador)
        k_medio = 0.5 * (indices[:, 2].min() + indices[:, 2].max())
        capa = indices[:, 2][np.argmin(np.abs(indices[:, 2] - k_medio))]
        en_capa = indices[indices[:, 2] == capa]
        eje_central = 0.5 * (indices[:, :2].min(axis=0) + indices[:, :2].max(axis=0))
        interior = tuple(int(i) for i in en_capa[
            np.argmin(np.sum((en_capa[:, :2] - eje_central) ** 2, axis=1))])
        for candidata in candidatas:
            eje, signo = "xyz".index(candidata[0]), (-1 if candidata[1] == "-" else 1)
            actual = list(interior)
            while True:
                siguiente = list(actual)
                siguiente[eje] += signo
                if (not 0 <= siguiente[eje] < grilla.forma[eje]
                        or grilla.vertida_id[tuple(siguiente)] != identificador):
                    break
                actual = siguiente
            if _vecino(grilla, tuple(actual), candidata) == 0:
                propuesta[str(identificador)] = {
                    "interior_m": [round(v, 6) for v in _centro_voxel(grilla, interior)],
                    "voxel_proximo_m": [round(v, 6) for v in _centro_voxel(grilla, tuple(actual))],
                    "cara_superficie": candidata,
                }
                break
        else:
            raise ValueError(f"Vertida {identificador}: no se encontró una cara lateral libre "
                             f"entre {candidatas} en la capa z={capa}")
    return propuesta


def temperatura_cara_estimada(condicion: CondicionBorde, temperatura_celda_c: float,
                              conductividad_w_m_k: float, ancho_celda_m: float,
                              tiempo_s: float, resistencia_extra_m2_k_w: float = 0.0) -> float:
    """Temperatura [°C] en la cara exterior de una celda con resistencia de media celda.

    ``resistencia_extra_m2_k_w`` es la protección (moldaje, aislación) en serie entre la
    cara del hormigón y el ambiente: la cara medida sigue siendo la del hormigón.
    """
    resistencia_celda = ancho_celda_m / (2.0 * conductividad_w_m_k)   # Rc = d/(2k) [m² K/W]
    if condicion.tipo == "dirichlet":
        return float(condicion.temperatura(tiempo_s))
    if condicion.tipo == "neumann":
        return float(temperatura_celda_c + condicion.flujo_entrante(tiempo_s) * resistencia_celda)
    h = condicion.coeficiente_h_w_m2_k
    if h <= 0.0:                     # Robin con h = 0 equivale a aislamiento
        return float(temperatura_celda_c)
    flujo_entrante = (condicion.temperatura(tiempo_s) - temperatura_celda_c) / (
        resistencia_celda + resistencia_extra_m2_k_w + 1.0 / h)
    return float(temperatura_celda_c + flujo_entrante * resistencia_celda)


def resumen_diferencias(tiempos_s, diferencia_c) -> dict | None:
    """Máximo |ΔT| sobre todas las muestras válidas, con signo e instante.

    Ignora NaN (n/a). Devuelve None si no hay ninguna muestra válida. En caso de
    empate conserva el primer instante.
    """
    tiempos = np.asarray(tiempos_s, dtype=float)
    diferencia = np.asarray(diferencia_c, dtype=float)
    validos = np.isfinite(diferencia)
    if not np.any(validos):
        return None
    absoluta = np.where(validos, np.abs(diferencia), -np.inf)
    posicion = int(np.argmax(absoluta))
    return {"max_abs_c": float(absoluta[posicion]),
            "diferencia_con_signo_c": float(diferencia[posicion]),
            "tiempo_h": float(tiempos[posicion] / 3600.0)}


class SeguimientoPuntos:
    """Registra T interior/superficial y estado de exposición en cada paso."""

    def __init__(self, grilla, puntos: list[PuntoVertida]):
        if not puntos:
            raise ValueError("Se requiere al menos un punto de seguimiento")
        self.grilla = grilla
        self.puntos = list(puntos)
        self.tiempos_s: list[float] = []
        self._datos = {p.vertida_id: {nombre: [] for nombre in (
            "interior", "celda", "cara", "estado")} for p in self.puntos}

    def registrar(self, tiempo_s: float, temperatura_c: np.ndarray, solver) -> None:
        """Toma una muestra. Si repite el último instante, lo reemplaza (activaciones)."""
        if self.tiempos_s and tiempo_s == self.tiempos_s[-1]:
            self.tiempos_s.pop()
            for datos in self._datos.values():
                for serie in datos.values():
                    serie.pop()
        self.tiempos_s.append(float(tiempo_s))
        for punto in self.puntos:
            datos = self._datos[punto.vertida_id]
            if not solver.mascara_activa[punto.indice_interior]:
                datos["interior"].append(math.nan)
                datos["celda"].append(math.nan)
                datos["cara"].append(math.nan)
                datos["estado"].append("n/a (vertida no activada)")
                continue
            t_celda = float(temperatura_c[punto.indice_voxel_proximo])
            datos["interior"].append(float(temperatura_c[punto.indice_interior]))
            datos["celda"].append(t_celda)
            expuesta = bool(solver.superficie_expuesta[punto.cara][punto.indice_voxel_proximo])
            vecino = _vecino(self.grilla, punto.indice_voxel_proximo, punto.cara)
            if expuesta:
                eje = "xyz".index(punto.cara[0])
                condicion, r_extra = solver.parametros_borde_celda(
                    punto.cara, punto.indice_voxel_proximo)
                datos["cara"].append(temperatura_cara_estimada(
                    condicion, t_celda,
                    float(solver.k[punto.indice_voxel_proximo]), solver.dx_m[eje], tiempo_s,
                    r_extra))
                datos["estado"].append("expuesta" if vecino in (0, punto.vertida_id)
                                       else f"expuesta (vecina futura: vertida {vecino})")
            else:
                datos["cara"].append(math.nan)
                datos["estado"].append(f"interfaz con vertida {vecino} (activa)")

    # -- resultados ---------------------------------------------------------
    def serie(self, vertida_id: int) -> dict:
        """Series numpy de una vertida (NaN = n/a) incluyendo ambas diferencias."""
        datos = self._datos[vertida_id]
        interior = np.asarray(datos["interior"], dtype=float)
        celda = np.asarray(datos["celda"], dtype=float)
        cara = np.asarray(datos["cara"], dtype=float)
        return {"tiempos_s": np.asarray(self.tiempos_s, dtype=float),
                "interior_c": interior, "celda_proxima_c": celda, "cara_estimada_c": cara,
                "dt_cara_c": interior - cara, "dt_celda_c": interior - celda,
                "estado": list(datos["estado"])}

    def resumen(self) -> list[dict]:
        """Una fila por vertida con el máximo |ΔT| (cara estimada y celda próxima)."""
        filas = []
        for punto in self.puntos:
            serie = self.serie(punto.vertida_id)
            filas.append({
                "vertida_id": punto.vertida_id,
                "punto": punto,
                "dt_cara": resumen_diferencias(serie["tiempos_s"], serie["dt_cara_c"]),
                "dt_celda": resumen_diferencias(serie["tiempos_s"], serie["dt_celda_c"]),
                "cambios_exposicion": self.cambios_de_estado(punto.vertida_id),
            })
        return filas

    def cambios_de_estado(self, vertida_id: int) -> list[tuple[float, str]]:
        """Instantes [h] del primer estado activo y de cada cambio posterior de exposición."""
        cambios, previo = [], None
        for tiempo, estado in zip(self.tiempos_s, self._datos[vertida_id]["estado"]):
            if estado.startswith("n/a"):
                continue
            if estado != previo:
                cambios.append((tiempo / 3600.0, estado))
                previo = estado
        return cambios

    def guardar_csv(self, ruta: Path) -> None:
        """CSV largo: una fila por instante de integración y vertida."""
        def fmt(valor):
            return NA if not math.isfinite(valor) else f"{valor:.8f}"
        with Path(ruta).open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow([
                "tiempo_h", "VertidaID",
                "x_centro_voxel_interior_m", "y_centro_voxel_interior_m", "z_centro_voxel_interior_m",
                "x_centro_voxel_proximo_m", "y_centro_voxel_proximo_m", "z_centro_voxel_proximo_m",
                "cara", "x_centro_cara_T_superficie_m", "y_centro_cara_T_superficie_m",
                "z_centro_cara_T_superficie_m",
                "T_interior_C", "T_celda_proxima_superficie_C", "T_cara_estimada_C",
                "dT_cara_C", "abs_dT_cara_C", "dT_celda_C", "abs_dT_celda_C",
                "estado_exposicion"])
            for punto in self.puntos:
                serie = self.serie(punto.vertida_id)
                fijos = [*(f"{v:.4f}" for v in punto.centro_interior_m),
                         *(f"{v:.4f}" for v in punto.centro_voxel_proximo_m),
                         punto.cara, *(f"{v:.4f}" for v in punto.centro_cara_m)]
                for n, tiempo in enumerate(serie["tiempos_s"]):
                    d_cara, d_celda = serie["dt_cara_c"][n], serie["dt_celda_c"][n]
                    escritor.writerow([
                        f"{tiempo / 3600.0:.6f}", punto.vertida_id, *fijos,
                        fmt(serie["interior_c"][n]), fmt(serie["celda_proxima_c"][n]),
                        fmt(serie["cara_estimada_c"][n]),
                        fmt(d_cara), fmt(abs(d_cara)) if math.isfinite(d_cara) else NA,
                        fmt(d_celda), fmt(abs(d_celda)) if math.isfinite(d_celda) else NA,
                        serie["estado"][n]])

    def guardar_resumen_csv(self, ruta: Path) -> None:
        def celda(resumen, clave, formato):
            return NA if resumen is None else format(resumen[clave], formato)
        with Path(ruta).open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow([
                "VertidaID",
                "x_centro_voxel_interior_m", "y_centro_voxel_interior_m", "z_centro_voxel_interior_m",
                "x_centro_voxel_proximo_m", "y_centro_voxel_proximo_m", "z_centro_voxel_proximo_m",
                "cara", "x_centro_cara_T_superficie_m", "y_centro_cara_T_superficie_m",
                "z_centro_cara_T_superficie_m", "distancia_centro_interior_a_cara_m",
                "max_abs_dT_cara_C", "dT_cara_con_signo_C", "hora_max_dT_cara_h",
                "max_abs_dT_celda_C", "dT_celda_con_signo_C", "hora_max_dT_celda_h",
                "cambios_de_exposicion"])
            for fila in self.resumen():
                punto = fila["punto"]
                cambios = "; ".join(f"{t:g} h: {estado}" for t, estado in fila["cambios_exposicion"])
                escritor.writerow([
                    punto.vertida_id,
                    *(f"{v:.4f}" for v in punto.centro_interior_m),
                    *(f"{v:.4f}" for v in punto.centro_voxel_proximo_m),
                    punto.cara, *(f"{v:.4f}" for v in punto.centro_cara_m),
                    f"{punto.distancia_interior_a_cara_m:.3f}",
                    celda(fila["dt_cara"], "max_abs_c", ".4f"),
                    celda(fila["dt_cara"], "diferencia_con_signo_c", "+.4f"),
                    celda(fila["dt_cara"], "tiempo_h", ".4f"),
                    celda(fila["dt_celda"], "max_abs_c", ".4f"),
                    celda(fila["dt_celda"], "diferencia_con_signo_c", "+.4f"),
                    celda(fila["dt_celda"], "tiempo_h", ".4f"),
                    cambios])
