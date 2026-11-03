"""Extracción de vertidas desde IFC4 mediante IfcOpenShell.

Ifcopenshell entrega coordenadas geométricas SI (metros) cuando
``convert-back-units`` está desactivado. La configuración también activa
``use-world-coords`` para incluir la colocación IFC de cada objeto.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np


@dataclass(frozen=True)
class VertidaIFC:
    """Malla triangulada y datos de planificación de una vertida."""

    vertida_id: int
    nombre: str
    global_id: str
    inicio_h: float
    temperatura_colocacion_c: float
    tipo_hormigon: str
    vertices_m: np.ndarray
    triangulos: np.ndarray
    propiedades: dict

    @property
    def minimo_m(self) -> np.ndarray:
        return self.vertices_m.min(axis=0)

    @property
    def maximo_m(self) -> np.ndarray:
        return self.vertices_m.max(axis=0)

    @property
    def dimensiones_m(self) -> np.ndarray:
        return self.maximo_m - self.minimo_m

    @property
    def volumen_m3(self) -> float:
        """Volumen firmado de la malla cerrada, expresado como magnitud."""
        a, b, c = (self.vertices_m[self.triangulos[:, i]] for i in range(3))
        return abs(float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)


def extraer_vertidas(ruta_ifc: str | Path) -> list[VertidaIFC]:
    """Lee las IfcFooting llamadas Vertida_XX y sus propiedades.

    No modifica el IFC. Falla claramente si falta geometría o una propiedad
    necesaria para identificar y activar una vertida.
    """
    ruta = Path(ruta_ifc).expanduser().resolve()
    if not ruta.is_file():
        raise FileNotFoundError(f"No existe el archivo IFC: {ruta}")

    modelo = ifcopenshell.open(str(ruta))
    configuracion = ifcopenshell.geom.settings()
    configuracion.set("use-world-coords", True)
    configuracion.set("convert-back-units", False)

    productos = [
        p for p in modelo.by_type("IfcProduct")
        if (getattr(p, "Name", "") or "").startswith("Vertida_")
    ]
    if not productos:
        raise ValueError("No se encontraron productos IFC con nombre Vertida_XX")

    resultado = []
    ids = set()
    for producto in productos:
        psets = ifcopenshell.util.element.get_psets(producto)
        props = dict(psets.get("UANDES_HormigonMasivo", {}))
        requeridas = ("VertidaID", "Inicio_h", "TemperaturaColocacion_C")
        faltantes = [nombre for nombre in requeridas if nombre not in props]
        if faltantes:
            raise ValueError(f"{producto.Name}: faltan propiedades IFC {faltantes}")

        vertida_id = int(props["VertidaID"])
        if vertida_id in ids:
            raise ValueError(f"VertidaID repetido: {vertida_id}")
        ids.add(vertida_id)

        forma = ifcopenshell.geom.create_shape(configuracion, producto).geometry
        vertices = np.asarray(forma.verts, dtype=float).reshape(-1, 3)
        triangulos = np.asarray(forma.faces, dtype=np.int64).reshape(-1, 3)
        if len(vertices) < 4 or len(triangulos) < 4:
            raise ValueError(f"{producto.Name}: la geometría no forma una malla cerrada útil")
        if not np.all(np.isfinite(vertices)):
            raise ValueError(f"{producto.Name}: la malla tiene coordenadas no finitas")

        resultado.append(VertidaIFC(
            vertida_id=vertida_id,
            nombre=producto.Name,
            global_id=producto.GlobalId,
            inicio_h=float(props["Inicio_h"]),
            temperatura_colocacion_c=float(props["TemperaturaColocacion_C"]),
            tipo_hormigon=str(props.get("TipoHormigon", "")),
            vertices_m=vertices,
            triangulos=triangulos,
            propiedades=props,
        ))

    return sorted(resultado, key=lambda vertida: vertida.vertida_id)
