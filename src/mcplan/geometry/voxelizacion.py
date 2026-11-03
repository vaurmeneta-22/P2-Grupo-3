"""Voxelización 3D de mallas IFC mediante prueba de punto interior."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ifc import VertidaIFC


@dataclass(frozen=True)
class GrillaVoxel:
    """Grilla regular. ``vertida_id=0`` identifica celdas fuera del sólido."""

    origen_m: tuple[float, float, float]
    paso_m: tuple[float, float, float]
    vertida_id: np.ndarray
    volumen_ifc_por_vertida_m3: dict[int, float]

    @property
    def forma(self) -> tuple[int, int, int]:
        return self.vertida_id.shape

    @property
    def volumen_voxel_m3(self) -> float:
        return float(np.prod(self.paso_m))

    @property
    def mascara_solido(self) -> np.ndarray:
        return self.vertida_id > 0

    @property
    def volumen_voxelizado_m3(self) -> float:
        return float(np.count_nonzero(self.mascara_solido) * self.volumen_voxel_m3)

    @property
    def volumen_voxelizado_por_vertida_m3(self) -> dict[int, float]:
        return {
            int(identificador): float(np.count_nonzero(self.vertida_id == identificador)
                                      * self.volumen_voxel_m3)
            for identificador in np.unique(self.vertida_id) if identificador > 0
        }

    @property
    def error_volumen_relativo_por_vertida(self) -> dict[int, float]:
        voxelizados = self.volumen_voxelizado_por_vertida_m3
        return {
            identificador: abs(voxelizados.get(identificador, 0.0) - volumen) / volumen
            for identificador, volumen in self.volumen_ifc_por_vertida_m3.items()
        }

    @property
    def error_volumen_relativo_total(self) -> float:
        volumen_ifc = sum(self.volumen_ifc_por_vertida_m3.values())
        if volumen_ifc <= 0:
            raise ValueError("El volumen total del IFC debe ser positivo")
        return abs(self.volumen_voxelizado_m3 - volumen_ifc) / volumen_ifc

    def coordenadas_centros_m(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ejes = [self.origen_m[i] + (np.arange(self.forma[i]) + 0.5) * self.paso_m[i]
                for i in range(3)]
        return tuple(ejes)


def _puntos_dentro_malla(puntos: np.ndarray, vertices: np.ndarray,
                         triangulos: np.ndarray) -> np.ndarray:
    """Clasifica puntos con paridad de intersecciones de rayo (malla cerrada)."""
    if len(puntos) == 0:
        return np.zeros(0, dtype=bool)

    direccion = np.array([1.0, 0.3713906763541037, 0.6947465906068658])
    direccion /= np.linalg.norm(direccion)
    a = vertices[triangulos[:, 0]]
    b = vertices[triangulos[:, 1]]
    c = vertices[triangulos[:, 2]]
    arista_1 = b - a
    arista_2 = c - a
    h = np.cross(np.broadcast_to(direccion, arista_2.shape), arista_2)
    determinante = np.einsum("ij,ij->i", arista_1, h)
    valida = np.abs(determinante) > 1e-12
    inverso = np.zeros_like(determinante)
    inverso[valida] = 1.0 / determinante[valida]

    intersecciones = np.zeros(len(puntos), dtype=np.int32)
    # Batches acotan memoria temporal para grillas finas.
    for inicio in range(0, len(puntos), 4096):
        lote = puntos[inicio:inicio + 4096]
        conteo = np.zeros(len(lote), dtype=np.int32)
        for indice in np.flatnonzero(valida):
            s = lote - a[indice]
            u = np.einsum("ij,j->i", s, h[indice]) * inverso[indice]
            dentro_u = (u >= -1e-10) & (u <= 1.0 + 1e-10)
            if not np.any(dentro_u):
                continue
            q = np.cross(s, arista_1[indice])
            v = (q @ direccion) * inverso[indice]
            t = (q @ arista_2[indice]) * inverso[indice]
            golpe = (dentro_u & (v >= -1e-10) & (u + v <= 1.0 + 1e-10)
                     & (t > 1e-10))
            conteo += golpe
        intersecciones[inicio:inicio + len(lote)] = conteo
    return (intersecciones % 2) == 1


def voxelizar_vertidas(vertidas: list[VertidaIFC], paso_m: float | tuple[float, float, float],
                       *, origen_m: tuple[float, float, float] | None = None) -> GrillaVoxel:
    """Asigna a cada centro de voxel el ID de la vertida que lo contiene.

    En caso de solape entre dos vertidas, lanza ValueError en lugar de ocultar
    la ambigüedad. La malla se considera cerrada y se usa una dirección oblicua
    de rayo para reducir aciertos degenerados en aristas.
    """
    if not vertidas:
        raise ValueError("Se requiere al menos una vertida")
    paso = ((float(paso_m),) * 3 if np.isscalar(paso_m)
            else tuple(float(d) for d in paso_m))
    if len(paso) != 3 or any(not np.isfinite(d) or d <= 0 for d in paso):
        raise ValueError("paso_m debe ser positivo y contener dx, dy y dz")

    minimo = np.min(np.stack([v.minimo_m for v in vertidas]), axis=0)
    maximo = np.max(np.stack([v.maximo_m for v in vertidas]), axis=0)
    origen = minimo if origen_m is None else np.asarray(origen_m, dtype=float)
    if origen.shape != (3,) or not np.all(np.isfinite(origen)) or np.any(origen > minimo + 1e-10):
        raise ValueError("origen_m debe ser un punto finito bajo el límite mínimo del modelo")
    forma = tuple(np.ceil((maximo - origen) / np.asarray(paso) - 1e-12).astype(int))
    ejes = [origen[i] + (np.arange(forma[i]) + 0.5) * paso[i] for i in range(3)]
    campo = np.zeros(forma, dtype=np.int32)
    vol_ifc = {}

    for vertida in vertidas:
        vol_ifc[vertida.vertida_id] = vertida.volumen_m3
        indices = []
        for eje in range(3):
            indices.append(np.flatnonzero((ejes[eje] >= vertida.minimo_m[eje] - 1e-10)
                                          & (ejes[eje] <= vertida.maximo_m[eje] + 1e-10)))
        if any(len(indice) == 0 for indice in indices):
            continue
        ix, iy, iz = np.meshgrid(*indices, indexing="ij")
        indexacion = (ix.ravel(), iy.ravel(), iz.ravel())
        puntos = np.column_stack([ejes[eje][indexacion[eje]] for eje in range(3)])
        interior = _puntos_dentro_malla(puntos, vertida.vertices_m, vertida.triangulos)
        celdas = tuple(componente[interior] for componente in indexacion)
        ocupadas = campo[celdas]
        solapadas = (ocupadas != 0) & (ocupadas != vertida.vertida_id)
        if np.any(solapadas):
            ids_solapados = np.unique(ocupadas[solapadas]).tolist()
            raise ValueError(f"Solape voxelizado en VertidaID={vertida.vertida_id} con {ids_solapados}")
        campo[celdas] = vertida.vertida_id

    return GrillaVoxel(tuple(float(x) for x in origen), paso, campo, vol_ifc)


def caras_superficie_exterior(vertida_id: np.ndarray,
                              mascara_activa: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """Máscaras 3D de caras expuestas a vacío o al borde de la grilla.

    Las máscaras tienen la forma de la grilla y True señala la cara de un
    voxel sólido que no comparte contacto con otro voxel activo.
    """
    ids = np.asarray(vertida_id)
    if ids.ndim != 3:
        raise ValueError("vertida_id debe ser un arreglo 3D")
    activos = ids > 0 if mascara_activa is None else np.asarray(mascara_activa, dtype=bool)
    if activos.shape != ids.shape:
        raise ValueError("mascara_activa debe tener la misma forma que vertida_id")
    expuestas = {}
    for eje, nombre in enumerate("xyz"):
        for lado in ("-", "+"):
            cara = np.zeros(ids.shape, dtype=bool)
            borde = [slice(None)] * 3
            borde[eje] = 0 if lado == "-" else -1
            cara[tuple(borde)] = activos[tuple(borde)]
            if ids.shape[eje] > 1:
                actual, vecino = [slice(None)] * 3, [slice(None)] * 3
                if lado == "-":
                    actual[eje], vecino[eje] = slice(1, None), slice(None, -1)
                else:
                    actual[eje], vecino[eje] = slice(None, -1), slice(1, None)
                cara[tuple(actual)] = activos[tuple(actual)] & ~activos[tuple(vecino)]
            expuestas[nombre + lado] = cara
    return expuestas


def area_interfaces_entre_vertidas(vertida_id: np.ndarray,
                                   paso_m: float | tuple[float, float, float],
                                   mascara_activa: np.ndarray | None = None) -> dict[tuple[int, int], float]:
    """Área de contacto voxelizada [m²] para cada par de vertidas vecinas."""
    ids = np.asarray(vertida_id)
    if ids.ndim != 3:
        raise ValueError("vertida_id debe ser un arreglo 3D")
    paso = ((float(paso_m),) * 3 if np.isscalar(paso_m)
            else tuple(float(d) for d in paso_m))
    if len(paso) != 3 or any(not np.isfinite(d) or d <= 0 for d in paso):
        raise ValueError("paso_m debe ser positivo y contener dx, dy y dz")
    activos = ids > 0 if mascara_activa is None else np.asarray(mascara_activa, dtype=bool)
    if activos.shape != ids.shape:
        raise ValueError("mascara_activa debe tener la misma forma que vertida_id")
    areas = {}
    for eje in range(3):
        bajo = [slice(None)] * 3
        alto = [slice(None)] * 3
        bajo[eje], alto[eje] = slice(None, -1), slice(1, None)
        bajo, alto = tuple(bajo), tuple(alto)
        ids_bajo, ids_alto = ids[bajo], ids[alto]
        contacto = (activos[bajo] & activos[alto] & (ids_bajo != ids_alto)
                    & (ids_bajo > 0) & (ids_alto > 0))
        pares, cuentas = np.unique(
            np.sort(np.stack((ids_bajo[contacto], ids_alto[contacto]), axis=1), axis=1),
            axis=0, return_counts=True)
        area = float(np.prod([paso[i] for i in range(3) if i != eje]))
        for par, cuenta in zip(pares, cuentas):
            clave = (int(par[0]), int(par[1]))
            areas[clave] = areas.get(clave, 0.0) + int(cuenta) * area
    return areas

