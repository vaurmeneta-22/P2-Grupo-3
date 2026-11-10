"""Exportación de estados 3D a VTK XML (.vtu + .pvd), legible en PyVista y ParaView.

No requiere VTK ni PyVista: escribe XML ASCII estándar.

Unidades y convenciones
-----------------------
* Coordenadas de los puntos: **metros** (mismo sistema que el IFC).
* Datos de celda (una celda = un voxel): ``temperatura_C`` [°C], ``vertida_id`` (entero),
  ``mascara_activa`` (1 = hormigón ya colocado), ``n_caras_expuestas`` (0 a 6: caras
  exteriores de ese voxel, sin vecino activo).
* Por defecto solo se escriben los voxeles ACTIVOS (``solo_activas=True``): el hormigón
  aún no colocado no existe en el archivo, así que la vista nunca lo muestra.
* Datos de campo (``FieldData``): ``TimeValue`` y ``tiempo_h`` [h], ``tiempo_s`` [s].
  El ``.pvd`` guarda el mismo tiempo en horas (``timestep``).
* Temperatura superficial: el dato de celda es la del CENTRO del voxel. La temperatura de la
  CARA exterior reconstruida (resistencia de media celda y condición de borde del solver) se
  exporta aparte, en ``marcas_*.vtu``, con ``codigo_marca = 3``. No se mezclan.

Marcas (``marcas_*.vtu``): puntos con ``codigo_marca``
    1 = núcleo (centro del voxel interior)
    2 = centro del voxel próximo a la cara exterior
    3 = centro de la cara exterior, temperatura RECONSTRUIDA
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

from mcplan.geometry.voxelizacion import caras_superficie_exterior

CODIGO_NUCLEO = 1
CODIGO_VOXEL_PROXIMO = 2
CODIGO_CARA_RECONSTRUIDA = 3
LEYENDA_MARCAS = {
    str(CODIGO_NUCLEO): "nucleo: centro del voxel interior [°C del voxel]",
    str(CODIGO_VOXEL_PROXIMO): "voxel_proximo_a_la_cara: centro del voxel de borde [°C del voxel]",
    str(CODIGO_CARA_RECONSTRUIDA): "cara_exterior_reconstruida: centro de la cara, T reconstruida [°C]",
}
VTK_HEXAHEDRON = 12
VTK_VERTEX = 1


@dataclass(frozen=True)
class Marca:
    vertida_id: int
    codigo: int
    xyz_m: tuple[float, float, float]
    temperatura_c: float


def _ascii(valores, por_linea=6) -> str:
    plano = np.asarray(valores).reshape(-1)
    if plano.dtype.kind == "f":
        piezas = [f"{v:.9g}" for v in plano]
    else:
        piezas = [str(int(v)) for v in plano]
    return "\n".join(" ".join(piezas[i:i + por_linea]) for i in range(0, len(piezas), por_linea))


def _array(nombre, tipo, valores, componentes=1, tuplas=None) -> str:
    comp = f' NumberOfComponents="{componentes}"' if componentes != 1 else ""
    if tuplas is not None:                      # FieldData exige NumberOfTuples para que VTK lea el dato
        comp += f' NumberOfTuples="{tuplas}"'
    return (f'<DataArray type="{tipo}" Name="{escape(nombre)}"{comp} format="ascii">\n'
            f"{_ascii(valores)}\n</DataArray>")


def _campo(tiempo_s: float) -> str:
    return ("<FieldData>\n"
            + _array("TimeValue", "Float64", [tiempo_s / 3600.0], tuplas=1) + "\n"
            + _array("tiempo_h", "Float64", [tiempo_s / 3600.0], tuplas=1) + "\n"
            + _array("tiempo_s", "Float64", [tiempo_s], tuplas=1) + "\n</FieldData>")


def _vtu(puntos, conectividad, offsets, tipos, datos_celda: str, datos_punto: str, tiempo_s) -> str:
    n_celdas = len(offsets)
    return f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">
<UnstructuredGrid>
{_campo(tiempo_s)}
<Piece NumberOfPoints="{len(puntos)}" NumberOfCells="{n_celdas}">
<PointData>
{datos_punto}
</PointData>
<CellData>
{datos_celda}
</CellData>
<Points>
{_array("coordenadas_m", "Float64", puntos, 3)}
</Points>
<Cells>
{_array("connectivity", "Int64", conectividad)}
{_array("offsets", "Int64", offsets)}
{_array("types", "UInt8", tipos)}
</Cells>
</Piece>
</UnstructuredGrid>
</VTKFile>
"""


def escribir_estado(ruta, grilla, temperatura_c, mascara_activa, tiempo_s, *, solo_activas=True) -> int:
    """Escribe un ``.vtu`` con un hexaedro por voxel. Devuelve el número de celdas."""
    ids = grilla.vertida_id
    seleccion = (ids > 0) & (np.asarray(mascara_activa, dtype=bool) if solo_activas else True)
    celdas = np.argwhere(seleccion)
    if len(celdas) == 0:
        raise ValueError("No hay voxeles para exportar en este estado")
    origen = np.asarray(grilla.origen_m, dtype=float)
    paso = np.asarray(grilla.paso_m, dtype=float)
    esquinas = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    indices = (celdas[:, None, :] + esquinas[None, :, :]).reshape(-1, 3)
    unicos, inverso = np.unique(indices, axis=0, return_inverse=True)
    puntos = origen + unicos * paso
    conectividad = np.asarray(inverso).reshape(-1)
    offsets = np.arange(1, len(celdas) + 1) * 8
    tipos = np.full(len(celdas), VTK_HEXAHEDRON)
    ii, jj, kk = celdas[:, 0], celdas[:, 1], celdas[:, 2]
    activa = np.asarray(mascara_activa, dtype=bool)
    caras = caras_superficie_exterior(ids, activa)
    n_expuestas = sum(c.astype(int) for c in caras.values())
    datos = "\n".join([
        _array("temperatura_C", "Float64", np.asarray(temperatura_c)[ii, jj, kk]),
        _array("vertida_id", "Int32", ids[ii, jj, kk]),
        _array("mascara_activa", "UInt8", activa[ii, jj, kk].astype(int)),
        _array("n_caras_expuestas", "UInt8", n_expuestas[ii, jj, kk]),
    ])
    Path(ruta).write_text(_vtu(puntos, conectividad, offsets, tipos, datos, "", tiempo_s),
                          encoding="utf-8")
    return len(celdas)


def escribir_marcas(ruta, marcas: list[Marca], tiempo_s) -> None:
    """Escribe los puntos marcados (núcleo, voxel próximo, cara reconstruida) como vértices."""
    if not marcas:
        raise ValueError("Se requiere al menos una marca")
    puntos = np.array([m.xyz_m for m in marcas], dtype=float)
    n = len(marcas)
    arreglos = "\n".join([
        _array("temperatura_C", "Float64", [m.temperatura_c for m in marcas]),
        _array("codigo_marca", "Int32", [m.codigo for m in marcas]),
        _array("vertida_id", "Int32", [m.vertida_id for m in marcas]),
    ])
    Path(ruta).write_text(
        _vtu(puntos, np.arange(n), np.arange(1, n + 1), np.full(n, VTK_VERTEX), arreglos, arreglos, tiempo_s),
        encoding="utf-8")


def escribir_pvd(ruta, entradas: list[tuple[float, str]]) -> None:
    """Colección temporal ``.pvd``: ``entradas`` = [(tiempo_h, archivo_relativo), ...]."""
    filas = "\n".join(f'    <DataSet timestep="{t:.6f}" group="" part="0" file="{escape(f)}"/>'
                      for t, f in entradas)
    Path(ruta).write_text('<?xml version="1.0"?>\n'
                          '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
                          f"  <Collection>\n{filas}\n  </Collection>\n</VTKFile>\n", encoding="utf-8")


def exportar_estados(carpeta, grilla, resultado, tiempos_h, *, marcas_por_estado=None,
                     prefijo="estado", solo_activas=True) -> dict:
    """Exporta los estados guardados más cercanos a ``tiempos_h``.

    ``marcas_por_estado``: función ``f(indice_estado, tiempo_s) -> list[Marca]`` (opcional).
    Devuelve un resumen con las rutas escritas. Nunca borra archivos previos: si ya existen
    se sobrescriben solo los de mismo nombre dentro de ``carpeta``.
    """
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    tiempos = np.asarray(resultado.tiempos_s, dtype=float)
    indices = sorted({int(np.argmin(np.abs(tiempos - float(h) * 3600.0))) for h in tiempos_h})
    estados, marcas_pvd, escritos = [], [], []
    for n in indices:
        t_s = float(tiempos[n])
        nombre = f"{prefijo}_t{t_s / 3600.0:07.2f}h.vtu"
        celdas = escribir_estado(carpeta / nombre, grilla, resultado.temperaturas_c[n],
                                 resultado.mascaras_activas[n], t_s, solo_activas=solo_activas)
        estados.append((t_s / 3600.0, nombre))
        escritos.append({"archivo": nombre, "tiempo_h": t_s / 3600.0, "celdas": celdas})
        if marcas_por_estado is not None:
            marcas = marcas_por_estado(n, t_s)
            if marcas:
                nombre_m = f"marcas_t{t_s / 3600.0:07.2f}h.vtu"
                escribir_marcas(carpeta / nombre_m, marcas, t_s)
                marcas_pvd.append((t_s / 3600.0, nombre_m))
    escribir_pvd(carpeta / f"{prefijo}.pvd", estados)
    if marcas_pvd:
        escribir_pvd(carpeta / "marcas.pvd", marcas_pvd)
    (carpeta / "leyenda_marcas.json").write_text(
        json.dumps({"codigo_marca": LEYENDA_MARCAS, "unidades": {
            "coordenadas": "m", "temperatura_C": "°C", "tiempo_h": "h"}}, indent=2, ensure_ascii=False),
        encoding="utf-8")
    return {"carpeta": str(carpeta), "pvd": str(carpeta / f"{prefijo}.pvd"),
            "pvd_marcas": str(carpeta / "marcas.pvd") if marcas_pvd else None, "estados": escritos}
