"""Loader meteorologico: entrega la temperatura ambiente variable que usa Robin.

Fuentes (``clima.fuente`` en el JSONC):

* ``constante``: un valor fijo (reproduce la configuracion anterior).
* ``sintetica``: serie determinista identificada como SINTETICA (demo y pruebas).
* ``snapshot``: serie normalizada de un snapshot guardado. Si no cubre el horizonte,
  informa el intervalo disponible y, solo si se pide ``extension``, completa la cola
  marcando esos puntos como ``extension:<metodo>`` (distinguibles del pronostico).

El tiempo del solver es ``t_s`` desde ``inicio_obra_utc`` (t = 0 de la simulacion).

Como incorporar viento y radiacion sin contar dos veces el mismo intercambio
(no implementado; ver docs/E10_entorno_y_clima.md):

* Viento: reemplaza el ``h`` constante por ``h(v) = h_natural + b v`` (un solo h
  convectivo-radiativo combinado por cara). No se suma un segundo coeficiente.
* Radiacion solar: se modela con una temperatura aire-sol
  ``T_as = T_aire + alfa G / h``, usada como temperatura ambiente del MISMO Robin; si en
  cambio se agrega un flujo Neumann solar, la temperatura de Robin queda en T_aire.
  Nunca ambas.
* Onda larga: si ``h`` ya es un coeficiente combinado (convectivo + radiativo
  linealizado), no se agrega un termino de onda larga aparte.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .serie import SerieTemperatura, a_utc, iso_utc
from .sintetica import AVISO_SINTETICA, serie_sintetica
from .sitio import Sitio, sitio_desde_config
from .snapshots import Cobertura, cobertura, extender_serie, leer_snapshot, listar_snapshots


@dataclass(frozen=True)
class AireVariable:
    """Temperatura ambiente como funcion del tiempo, con su procedencia."""

    funcion: Callable[[float], float]
    descripcion: dict = field(default_factory=dict)
    serie: SerieTemperatura | None = None
    cobertura: Cobertura | None = None

    @property
    def es_sintetica(self) -> bool:
        return bool(self.descripcion.get("sintetica", False))


def construir_aire(config_clima: dict | None, *, tiempo_final_s: float, raiz_repo: Path,
                   temperatura_constante_c: float = 20.0) -> AireVariable:
    """Crea el ambiente segun ``config_clima``; sin bloque, usa temperatura constante."""
    cfg = config_clima or {}
    fuente = str(cfg.get("fuente", "constante")).lower()
    if fuente == "constante":
        valor = float((cfg.get("constante") or {}).get("temperatura_c", temperatura_constante_c))
        return AireVariable(lambda _t, v=valor: v,
                            {"fuente": "constante", "temperatura_c": valor, "sintetica": False})

    if fuente not in ("sintetica", "snapshot"):
        raise ValueError("clima.fuente debe ser 'constante', 'sintetica' o 'snapshot'")
    if not cfg.get("inicio_obra_utc") and fuente == "sintetica":
        raise ValueError("clima.inicio_obra_utc es obligatorio para la fuente sintetica")
    horas_necesarias = tiempo_final_s / 3600.0

    if fuente == "sintetica":
        inicio = a_utc(cfg["inicio_obra_utc"])
        parametros = cfg.get("sintetica") or {}
        serie = serie_sintetica(inicio, horas_necesarias + 1.0, **{
            k: float(v) for k, v in parametros.items()
            if k in ("media_c", "amplitud_diaria_c", "hora_pico_utc", "tendencia_c_por_dia", "paso_h")})
        descripcion = {"fuente": "sintetica", "sintetica": True, "aviso": AVISO_SINTETICA,
                       "parametros": parametros, "inicio_obra_utc": iso_utc(inicio)}
        return AireVariable(serie.a_funcion(inicio), descripcion, serie,
                            cobertura(serie, inicio, inicio_mas(inicio, tiempo_final_s)))

    cfg_snap = cfg.get("snapshot") or {}
    carpeta = cfg_snap.get("carpeta")
    if carpeta is None:
        existentes = listar_snapshots(Path(raiz_repo) / cfg_snap.get("raiz", "data/weather"))
        if not existentes:
            raise FileNotFoundError("clima.snapshot: no hay snapshots guardados en data/weather")
        carpeta = existentes[-1]
    else:
        carpeta = Path(raiz_repo) / carpeta
    instantanea = leer_snapshot(carpeta)
    # Sin inicio_obra_utc, la obra comienza en el primer instante cubierto por el snapshot.
    inicio = a_utc(cfg["inicio_obra_utc"]) if cfg.get("inicio_obra_utc") else instantanea.serie.intervalo[0]
    fin = inicio_mas(inicio, tiempo_final_s)
    cob = cobertura(instantanea.serie, inicio, fin)
    serie = instantanea.serie
    extension = cfg_snap.get("extension")
    descripcion = {"fuente": "snapshot", "sintetica": bool(instantanea.metadatos.get("synthetic")),
                   "tipo": instantanea.metadatos.get("series_kind", "pronostico"),
                   "snapshot": carpeta.name, "modelo": instantanea.metadatos.get("model"),
                   "descargado_utc": instantanea.metadatos.get("retrieved_at"),
                   "cobertura": cob.resumen(), "extension": extension,
                   "inicio_obra_utc": iso_utc(inicio)}
    if cob.horas_faltantes_inicio > 0:
        raise ValueError(f"El snapshot no cubre el inicio de obra: {cob.resumen()}")
    if not cob.cubre_completo:
        if not extension:
            raise ValueError("El snapshot no cubre el horizonte y no se pidio extension. "
                             + cob.resumen())
        serie = extender_serie(serie, fin, str(extension))
        descripcion["puntos_por_origen"] = serie.fraccion_por_origen()
        descripcion["extension_desde_utc"] = iso_utc(instantanea.serie.intervalo[1])
    return AireVariable(serie.a_funcion(inicio), descripcion, serie, cob)


def inicio_mas(inicio, segundos: float):
    from datetime import timedelta
    return a_utc(inicio) + timedelta(seconds=float(segundos))
