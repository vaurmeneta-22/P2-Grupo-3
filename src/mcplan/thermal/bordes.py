"""Condiciones de borde del entorno: aire, suelo y protecciones (moldaje, aislación).

Qué distingue este módulo (E10)
-------------------------------
* **Aire**: toda cara exterior expuesta usa Robin con ``h`` configurable y una
  temperatura ambiente que puede variar en el tiempo (ver ``mcplan.weather``).
* **Suelo**: la base física (celdas cuya cara ``z-`` está en el plano más bajo del
  hormigón) usa un modelo EFECTIVO: una capa de suelo de espesor ``L`` y
  conductividad ``k_s`` hacia una temperatura profunda ``T_suelo``. Equivale a un
  Robin con ``h_eq = k_s / L``. Supuestos (documentados en docs/E10_entorno_y_clima.md):
  régimen casi estacionario en el suelo (sin acumulación de calor en él), temperatura
  profunda conocida y uniforme, sin humedad ni napa. Un voladizo con cara ``z-``
  expuesta que no esté en la base usa aire.
* **Protecciones** (moldaje, aislación): capas con resistencia ``R'' = e / k``
  [m² K/W] que se suman EN SERIE con la media celda de hormigón y la convección:

      R_total = d/(2 k_hormigón) + Σ R''_capas + 1/h        (flujo = ΔT / R_total)

  Cada protección se define por vertida y cara, y puede retirarse. El retiro se
  mide explícitamente desde ``colocacion`` (hora de inicio de ESA vertida, por
  defecto) o desde ``inicio_obra`` (t = 0 de la simulación).

El intercambio por conducción entre celdas activas no se toca: estas condiciones
solo actúan en caras expuestas (sin vecino activo).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .diffusion_nd import CondicionBorde, robin

CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")
CARAS_LATERALES = ("x-", "x+", "y-", "y+")
REFERENCIAS_RETIRO = ("colocacion", "inicio_obra")


@dataclass(frozen=True)
class ModeloSuelo:
    """Suelo efectivo bajo la base: Robin con ``h_eq = k_s / L`` hacia ``T_suelo``."""

    temperatura_c: float | Callable[[float], float]
    conductividad_w_m_k: float
    espesor_efectivo_m: float

    def __post_init__(self):
        if not (self.conductividad_w_m_k > 0 and self.espesor_efectivo_m > 0):
            raise ValueError("suelo: conductividad y espesor efectivo deben ser positivos")

    @property
    def resistencia_m2_k_w(self) -> float:
        return self.espesor_efectivo_m / self.conductividad_w_m_k

    @property
    def h_equivalente_w_m2_k(self) -> float:
        return 1.0 / self.resistencia_m2_k_w

    def condicion(self) -> CondicionBorde:
        return robin(self.h_equivalente_w_m2_k, self.temperatura_c)


@dataclass(frozen=True)
class Proteccion:
    """Capa de moldaje o aislación sobre caras de una vertida."""

    vertida_id: int
    nombre: str
    caras: tuple[str, ...]
    resistencia_m2_k_w: float
    retiro_s: float | None          # instante absoluto (desde t=0 de la simulación); None = nunca
    referencia_retiro: str = "colocacion"

    def activa_en(self, tiempo_s: float) -> bool:
        return self.retiro_s is None or tiempo_s < self.retiro_s - 1e-9


def resistencia_de_capa(cfg: dict, contexto: str = "protección") -> float:
    """R'' [m² K/W] desde ``resistencia_m2_k_w`` o desde ``espesor_m`` y ``conductividad_w_m_k``."""
    if cfg.get("resistencia_m2_k_w") is not None:
        resistencia = float(cfg["resistencia_m2_k_w"])
    elif cfg.get("espesor_m") is not None and cfg.get("conductividad_w_m_k") is not None:
        conductividad = float(cfg["conductividad_w_m_k"])
        if conductividad <= 0:
            raise ValueError(f"{contexto}: conductividad_w_m_k debe ser positiva")
        resistencia = float(cfg["espesor_m"]) / conductividad
    else:
        raise ValueError(f"{contexto}: define resistencia_m2_k_w o espesor_m + conductividad_w_m_k")
    if not np.isfinite(resistencia) or resistencia < 0:
        raise ValueError(f"{contexto}: la resistencia debe ser finita y no negativa")
    return resistencia


def proteccion_desde_config(vertida, cfg: dict) -> Proteccion:
    """Convierte una entrada de configuración en una :class:`Proteccion` validada.

    ``retiro_h`` puede ser un número [h], ``null`` (no se retira) o ``"ifc"`` (usa la
    propiedad IFC ``RetiroMoldaje_h``; su referencia temporal no está documentada en el
    IFC, por eso se interpreta con ``retiro_desde``, por defecto ``colocacion``).
    """
    nombre = str(cfg.get("nombre", "proteccion"))
    contexto = f"Vertida {vertida.vertida_id} · {nombre}"
    caras = tuple(cfg.get("caras", CARAS_LATERALES))
    if not caras or any(c not in CARAS for c in caras):
        raise ValueError(f"{contexto}: caras debe ser una lista dentro de {CARAS}")
    referencia = str(cfg.get("retiro_desde", "colocacion")).lower()
    if referencia not in REFERENCIAS_RETIRO:
        raise ValueError(f"{contexto}: retiro_desde debe ser uno de {REFERENCIAS_RETIRO}")
    retiro = cfg.get("retiro_h")
    if isinstance(retiro, str):
        if retiro.lower() != "ifc":
            raise ValueError(f"{contexto}: retiro_h debe ser número, null o 'ifc'")
        if "RetiroMoldaje_h" not in vertida.propiedades:
            raise ValueError(f"{contexto}: el IFC no tiene RetiroMoldaje_h")
        retiro = float(vertida.propiedades["RetiroMoldaje_h"])
    if retiro is None:
        retiro_s = None
    else:
        retiro = float(retiro)
        if retiro < 0:
            raise ValueError(f"{contexto}: retiro_h no puede ser negativo")
        base_h = vertida.inicio_h if referencia == "colocacion" else 0.0
        retiro_s = (base_h + retiro) * 3600.0
        if retiro_s < vertida.inicio_h * 3600.0:
            raise ValueError(f"{contexto}: el retiro ocurre antes de la colocación de la vertida")
    return Proteccion(int(vertida.vertida_id), nombre, caras,
                      resistencia_de_capa(cfg, contexto), retiro_s, referencia)


class ModeloBordes:
    """Reúne aire, suelo y protecciones y las entrega al solver por cara y celda."""

    def __init__(self, grilla, condicion_aire: CondicionBorde, *, suelo: ModeloSuelo | None = None,
                 protecciones: tuple[Proteccion, ...] = ()):
        self.grilla = grilla
        self.condicion_aire = condicion_aire
        self.suelo = suelo
        self.protecciones = tuple(protecciones)
        ids_presentes = {int(i) for i in np.unique(grilla.vertida_id) if i > 0}
        for proteccion in self.protecciones:
            if proteccion.vertida_id not in ids_presentes:
                raise ValueError(f"Protección '{proteccion.nombre}': la vertida "
                                 f"{proteccion.vertida_id} no está en la grilla")
        solido = grilla.vertida_id > 0
        capas = np.flatnonzero(solido.any(axis=(0, 1)))
        self.indice_base = int(capas[0]) if capas.size else 0
        self.mascara_base = solido & (np.arange(grilla.forma[2])[None, None, :] == self.indice_base)

    # -- lo que consume el solver ------------------------------------------------
    def fronteras(self) -> dict[str, CondicionBorde]:
        """Condición por defecto de cada cara (aire)."""
        return {cara: self.condicion_aire for cara in CARAS}

    def fronteras_celda(self) -> dict:
        """Sobreescrituras por celda: la base física en contacto con el suelo."""
        if self.suelo is None:
            return {}
        return {"z-": [(self.mascara_base, self.suelo.condicion())]}

    def resistencias_en(self, tiempo_s: float) -> dict[str, np.ndarray]:
        """R'' extra [m² K/W] por cara y celda con las protecciones vigentes en ``tiempo_s``."""
        resultado: dict[str, np.ndarray] = {}
        for proteccion in self.protecciones:
            if not proteccion.activa_en(tiempo_s):
                continue
            celdas = self.grilla.vertida_id == proteccion.vertida_id
            for cara in proteccion.caras:
                campo = resultado.setdefault(cara, np.zeros(self.grilla.forma, dtype=float))
                campo[celdas] += proteccion.resistencia_m2_k_w     # capas en serie: se suman
        return resultado

    def instantes_retiro_s(self, tiempo_final_s: float) -> list[float]:
        """Instantes (0 < t < t_final) en que se retira al menos una protección."""
        return sorted({p.retiro_s for p in self.protecciones
                       if p.retiro_s is not None and 0.0 < p.retiro_s < tiempo_final_s})

    def retiradas_en(self, tiempo_s: float) -> list[Proteccion]:
        return [p for p in self.protecciones
                if p.retiro_s is not None and abs(p.retiro_s - tiempo_s) <= 1e-6]

    def aplicar(self, solver, tiempo_s: float) -> None:
        """Actualiza el solver con las protecciones vigentes; no modifica la temperatura."""
        solver.actualizar_bordes(resistencias_borde=self.resistencias_en(tiempo_s),
                                 fronteras_celda=self.fronteras_celda())


def crear_modelo_bordes(grilla, vertidas, config_entorno: dict | None,
                        condicion_aire: CondicionBorde) -> ModeloBordes:
    """Construye el modelo desde el bloque ``entorno`` del JSONC.

    ``entorno.suelo``: ``activo``, ``temperatura_c``, ``conductividad_w_m_k``,
    ``espesor_efectivo_m``. ``entorno.protecciones_por_vertida``: ``{"1": [ {...}, ... ]}``
    con ``nombre``, ``caras``, ``espesor_m`` + ``conductividad_w_m_k`` (o
    ``resistencia_m2_k_w``), ``retiro_h`` y ``retiro_desde``.
    """
    cfg = config_entorno or {}
    suelo = None
    cfg_suelo = cfg.get("suelo") or {}
    if cfg_suelo.get("activo", False):
        suelo = ModeloSuelo(
            temperatura_c=float(cfg_suelo.get("temperatura_c", 18.0)),
            conductividad_w_m_k=float(cfg_suelo.get("conductividad_w_m_k", 1.5)),
            espesor_efectivo_m=float(cfg_suelo.get("espesor_efectivo_m", 0.8)))
    por_id = {int(v.vertida_id): v for v in vertidas}
    protecciones = []
    for clave, lista in sorted((cfg.get("protecciones_por_vertida") or {}).items(),
                               key=lambda par: int(par[0])):
        if int(clave) not in por_id:
            raise ValueError(f"protecciones_por_vertida: no existe la vertida {clave}")
        for entrada in lista or []:
            protecciones.append(proteccion_desde_config(por_id[int(clave)], entrada))
    return ModeloBordes(grilla, condicion_aire, suelo=suelo, protecciones=tuple(protecciones))
