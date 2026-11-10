"""Sitio climático asignado al grupo (coordenadas exactas pendientes de sitios.yaml)."""
from __future__ import annotations

from dataclasses import dataclass


class SitioSinCoordenadas(RuntimeError):
    """No se puede consultar un servicio real sin latitud y longitud verificadas."""


@dataclass(frozen=True)
class Sitio:
    grupo: int
    nombre: str
    region: str
    latitud: float | None
    longitud: float | None
    fuente_coordenadas: str

    @property
    def tiene_coordenadas(self) -> bool:
        return self.latitud is not None and self.longitud is not None

    def exigir_coordenadas(self) -> tuple[float, float]:
        if not self.tiene_coordenadas:
            raise SitioSinCoordenadas(
                f"Sitio del Grupo {self.grupo} ({self.nombre}, {self.region}): faltan latitud y "
                f"longitud. {self.fuente_coordenadas}. Completa 'sitio' en la configuración con "
                "las coordenadas exactas de sitios.yaml antes de consultar Open-Meteo.")
        lat, lon = float(self.latitud), float(self.longitud)
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError("latitud/longitud fuera de rango")
        return lat, lon


def sitio_desde_config(cfg: dict | None) -> Sitio:
    cfg = cfg or {}
    return Sitio(
        grupo=int(cfg.get("grupo", 3)),
        nombre=str(cfg.get("nombre", "Ovalle / Ecos del Mar")),
        region=str(cfg.get("region", "Coquimbo")),
        latitud=None if cfg.get("latitud") is None else float(cfg["latitud"]),
        longitud=None if cfg.get("longitud") is None else float(cfg["longitud"]),
        fuente_coordenadas=str(cfg.get("fuente_coordenadas",
                                       "Coordenadas pendientes de sitios.yaml")))
