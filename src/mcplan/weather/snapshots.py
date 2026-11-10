"""Snapshots de pronóstico: archivo original, serie normalizada y metadatos de procedencia.

Cada descarga crea una carpeta NUEVA ``<raiz>/<descarga UTC>/`` (nunca se sobrescribe
ni se borra una anterior; si el nombre existe se agrega ``_2``, ``_3``, ...):

    original.<ext>   bytes exactos recibidos (sha256 en metadata.json)
    forecast.json    serie normalizada: time_utc (ISO, UTC), temperature_2m [degC]
    metadata.json    fuente, modelo, coordenadas, emisión (si existe), descarga,
                     validez, unidades, sha256 del original y ``series_kind``
                     (pronostico | reanalisis | observacion | sintetica)

La cobertura se informa aparte (``cobertura``) y toda extensión del pronóstico original
(``extender_serie``) marca cada punto agregado como ``extension:<método>``.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path

import numpy as np

from .serie import (ORIGEN_PRONOSTICO, ORIGEN_SINTETICA, TIPOS_SERIE, SerieTemperatura, a_utc,
                    iso_utc)

UNIDADES = {"temperature_2m": "degC"}


@dataclass(frozen=True)
class Snapshot:
    carpeta: Path
    serie: SerieTemperatura
    metadatos: dict


@dataclass(frozen=True)
class Cobertura:
    solicitado: tuple[datetime, datetime]
    disponible: tuple[datetime, datetime]
    cubre_completo: bool
    horas_faltantes_inicio: float
    horas_faltantes_final: float

    def resumen(self) -> str:
        d0, d1 = (iso_utc(t) for t in self.disponible)
        s0, s1 = (iso_utc(t) for t in self.solicitado)
        if self.cubre_completo:
            return f"cobertura completa: disponible {d0} a {d1} cubre {s0} a {s1}"
        return (f"cobertura INCOMPLETA: disponible {d0} a {d1}; solicitado {s0} a {s1}; "
                f"faltan {self.horas_faltantes_inicio:.1f} h al inicio y "
                f"{self.horas_faltantes_final:.1f} h al final")


def _carpeta_nueva(raiz: Path, descargado_utc: datetime) -> Path:
    base = a_utc(descargado_utc).strftime("%Y-%m-%dT%H%M%SZ")
    raiz.mkdir(parents=True, exist_ok=True)
    for intento in range(1, 1000):
        carpeta = raiz / (base if intento == 1 else f"{base}_{intento}")
        try:
            carpeta.mkdir()                 # falla si existe: conserva descargas anteriores
            return carpeta
        except FileExistsError:
            continue
    raise RuntimeError("no se pudo crear una carpeta de snapshot nueva")


def guardar_snapshot(raiz, *, bruto: bytes, nombre_original: str, serie: SerieTemperatura,
                     fuente: str, modelo: str | None, latitud: float | None,
                     longitud: float | None, zona_horaria: str = "UTC",
                     descargado_utc: datetime, emitido_utc: datetime | None = None,
                     sitio: str | None = None, tipo: str = "pronostico",
                     procedencia_coordenadas: str | None = None,
                     extra: dict | None = None) -> Path:
    """Guarda un snapshot nuevo y devuelve su carpeta."""
    if tipo not in TIPOS_SERIE:
        raise ValueError(f"tipo debe ser uno de {TIPOS_SERIE}")
    carpeta = _carpeta_nueva(Path(raiz), descargado_utc)
    (carpeta / nombre_original).write_bytes(bruto)
    inicio, fin = serie.intervalo
    forecast = {"time_utc": [iso_utc(t) for t in serie.tiempos_utc],
                "temperature_2m": [round(float(v), 4) for v in serie.temperatura_c],
                "units": UNIDADES}
    (carpeta / "forecast.json").write_text(json.dumps(forecast, indent=2), encoding="utf-8")
    metadatos = {
        "source": fuente,
        "model": modelo,
        "site": sitio,
        "latitude": latitud,
        "longitude": longitud,
        "timezone": zona_horaria,
        "coordinates_source": procedencia_coordenadas,
        "issued_at": None if emitido_utc is None else iso_utc(emitido_utc),
        "issued_at_note": ("emision del modelo no informada por la fuente"
                           if emitido_utc is None else "emision informada por la fuente"),
        "retrieved_at": iso_utc(descargado_utc),
        "valid_from": iso_utc(inicio),
        "valid_to": iso_utc(fin),
        "variables": ["temperature_2m"],
        "units": UNIDADES,
        "original_file": nombre_original,
        "original_sha256": hashlib.sha256(bruto).hexdigest(),
        "n_points": len(serie.tiempos_utc),
        "series_kind": tipo,          # pronostico | reanalisis | observacion | sintetica
        "synthetic": tipo == "sintetica",
    }
    if extra:
        metadatos["extra"] = extra
    (carpeta / "metadata.json").write_text(json.dumps(metadatos, indent=2, ensure_ascii=False),
                                           encoding="utf-8")
    return carpeta


def leer_snapshot(carpeta) -> Snapshot:
    """Lee y valida un snapshot (unidades, integridad del original y orden temporal)."""
    carpeta = Path(carpeta)
    metadatos = json.loads((carpeta / "metadata.json").read_text(encoding="utf-8"))
    forecast = json.loads((carpeta / "forecast.json").read_text(encoding="utf-8"))
    if metadatos.get("units", {}).get("temperature_2m") != "degC":
        raise ValueError(f"{carpeta.name}: la temperatura debe estar en degC")
    original = carpeta / metadatos["original_file"]
    if hashlib.sha256(original.read_bytes()).hexdigest() != metadatos["original_sha256"]:
        raise ValueError(f"{carpeta.name}: el archivo original no coincide con su sha256")
    origen = metadatos.get("series_kind") or (
        ORIGEN_SINTETICA if metadatos.get("synthetic") else ORIGEN_PRONOSTICO)
    tiempos = tuple(a_utc(t) for t in forecast["time_utc"])
    serie = SerieTemperatura(tiempos, np.asarray(forecast["temperature_2m"], dtype=float),
                             tuple(origen for _ in tiempos))
    return Snapshot(carpeta, serie, metadatos)


def listar_snapshots(raiz) -> list[Path]:
    """Carpetas de snapshot existentes, ordenadas por nombre (fecha de descarga)."""
    raiz = Path(raiz)
    if not raiz.is_dir():
        return []
    return sorted(c for c in raiz.iterdir() if (c / "metadata.json").is_file())


def cobertura(serie: SerieTemperatura, inicio_utc, fin_utc) -> Cobertura:
    """Compara el intervalo disponible con el solicitado (no inventa ni rellena datos)."""
    s0, s1 = a_utc(inicio_utc), a_utc(fin_utc)
    if s1 <= s0:
        raise ValueError("el intervalo solicitado debe tener duracion positiva")
    d0, d1 = serie.intervalo
    falta_ini = max((d0 - s0).total_seconds(), 0.0) / 3600.0
    falta_fin = max((s1 - d1).total_seconds(), 0.0) / 3600.0
    return Cobertura((s0, s1), (d0, d1), falta_ini == 0.0 and falta_fin == 0.0,
                     falta_ini, falta_fin)


def extender_serie(serie: SerieTemperatura, hasta_utc, metodo: str = "ciclo_diario") -> SerieTemperatura:
    """Completa la cola faltante hasta ``hasta_utc``, marcando cada punto como extension.

    * ``ciclo_diario``: repite el ultimo ciclo de 24 h del pronostico (requiere >= 24 h).
    * ``ultimo_valor``: mantiene constante la ultima temperatura.

    Es una hipotesis del equipo, NO un pronostico: por eso queda registrada en ``origen``.
    """
    if metodo not in ("ciclo_diario", "ultimo_valor"):
        raise ValueError("metodo debe ser 'ciclo_diario' o 'ultimo_valor'")
    hasta = a_utc(hasta_utc)
    fin = serie.tiempos_utc[-1]
    if hasta <= fin:
        return serie
    if metodo == "ciclo_diario" and (fin - serie.tiempos_utc[0]) < timedelta(hours=24):
        raise ValueError("ciclo_diario requiere al menos 24 h de pronostico")
    paso = serie.tiempos_utc[-1] - serie.tiempos_utc[-2]
    tiempos = list(serie.tiempos_utc)
    valores = [float(v) for v in serie.temperatura_c]
    origen = list(serie.origen)
    etiqueta = f"extension:{metodo}"
    t = fin + paso
    while t <= hasta + timedelta(seconds=1):
        if metodo == "ciclo_diario":
            # Mismo instante del dia, tomado del ultimo dia disponible del pronostico.
            dias = int(np.ceil((t - fin) / timedelta(days=1)))
            referencia = t - timedelta(days=dias)
            previos = [k for k in serie.tiempos_utc if k <= referencia]
            indice = serie.tiempos_utc.index(max(previos)) if previos else len(serie.tiempos_utc) - 1
            valor = float(serie.temperatura_c[indice])
        else:
            valor = float(serie.temperatura_c[-1])
        tiempos.append(t)
        valores.append(valor)
        origen.append(etiqueta)
        t += paso
    return SerieTemperatura(tuple(tiempos), np.asarray(valores), tuple(origen))
