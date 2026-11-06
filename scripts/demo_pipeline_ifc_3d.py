"""Pipeline IFC -> voxeles -> activación por vertidas -> hidratación 3D.

Edición rápida: cambia configs/semana3_fundacion.jsonc. El IFC se abre solo en
lectura y la secuencia se toma de sus propiedades salvo overrides en el JSON.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from mcplan.geometry.ifc import VertidaIFC, extraer_vertidas
from mcplan.geometry.voxelizacion import (
    area_interfaces_entre_vertidas,
    caras_superficie_exterior,
    voxelizar_vertidas,
)
from mcplan.thermal.analisis import (
    PuntoAnalisis,
    cortes_por_punto,
    eventos_activacion,
    metricas_tmax,
    seleccionar_puntos_representativos,
    serie_punto,
)
from mcplan.thermal.diffusion_nd import aislado, dirichlet, neumann, robin
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.vertidas import ResultadoEtapas, simular_vertidas


CONFIG_PREDETERMINADA = RAIZ / "configs" / "semana3_fundacion.jsonc"
SALIDA = RAIZ / "results" / "semana3"
SALIDA_E09 = RAIZ / "results" / "semana3_e09"
EVIDENCIA_E09 = RAIZ / "reports" / "semana3" / "evidencia_E09"
COLORES_VERTIDA = {1: "#2c7fb8", 2: "#41b6c4", 3: "#a1dab4", 4: "#fed976", 5: "#f03b20"}
COLORES_VERTIDA_VISTAS = {
    1: "#2c7fb8",
    2: "#2aa7b5",
    3: "#4daf4a",
    4: "#e8a33d",
    5: "#d7301f",
}


def argumentos():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG_PREDETERMINADA,
                        help="JSONC de entradas (default: configs/semana3_fundacion.jsonc)")
    parser.add_argument("--dx", type=float, help="Sobrescribe el paso de voxel en metros")
    parser.add_argument("--horas", type=float, help="Sobrescribe el horizonte de simulación")
    parser.add_argument("--mostrar-vertidas", type=str,
                        help="IDs separados por coma para la imagen geométrica; ejemplo: 1,2")
    parser.add_argument("--salida", type=Path, default=None,
                        help="Carpeta de resultados (default: results/semana3; "
                             "con --modo-e09: results/semana3_e09)")
    parser.add_argument("--evidencia", type=Path, default=None,
                        help="Carpeta versionable de evidencia E09 "
                             "(con --modo-e09: misma carpeta que --salida)")
    parser.add_argument("--modo-e09", action="store_true",
                        help="Genera solo las cinco figuras E09 y sus CSV en "
                             "carpetas exclusivas; omite GIF y figuras generales")
    parser.add_argument("--sin-gif", action="store_true",
                        help="Omite los GIF para acelerar una corrida de diagnóstico")
    parser.add_argument("--mostrar", action="store_true", help="Abre las figuras al terminar")
    return parser.parse_args()


def archivar_contenido(carpeta: Path, destino: Path) -> list[Path]:
    """Mueve los archivos de ``carpeta`` a ``destino`` sin borrar nada.

    Se usa en el modo E09 para que la carpeta exclusiva de la corrida quede
    limpia antes de generar las evidencias nuevas. No toca subdirectorios ni
    carpetas de otras semanas.
    """
    if not carpeta.exists():
        return []
    archivos = [ruta for ruta in carpeta.iterdir() if ruta.is_file()]
    if not archivos:
        return []
    destino.mkdir(parents=True, exist_ok=True)
    movidos = []
    for ruta in archivos:
        objetivo = destino / ruta.name
        if objetivo.exists():
            objetivo = destino / f"{ruta.stem}_{int(ruta.stat().st_mtime)}{ruta.suffix}"
        ruta.replace(objetivo)
        movidos.append(objetivo)
    return movidos


def quitar_comentarios_jsonc(texto: str) -> str:
    """Elimina comentarios // fuera de cadenas para leer JSONC con json.loads."""
    salida = []
    en_cadena = False
    escapado = False
    indice = 0
    while indice < len(texto):
        caracter = texto[indice]
        if en_cadena:
            salida.append(caracter)
            if escapado:
                escapado = False
            elif caracter == "\\":
                escapado = True
            elif caracter == '"':
                en_cadena = False
            indice += 1
        elif caracter == '"':
            en_cadena = True
            salida.append(caracter)
            indice += 1
        elif caracter == "/" and indice + 1 < len(texto) and texto[indice + 1] == "/":
            while indice < len(texto) and texto[indice] not in "\r\n":
                indice += 1
        else:
            salida.append(caracter)
            indice += 1
    return "".join(salida)


def cargar_configuracion(ruta: Path) -> dict:
    ruta = ruta if ruta.is_absolute() else RAIZ / ruta
    with ruta.open("r", encoding="utf-8") as archivo:
        config = json.loads(quitar_comentarios_jsonc(archivo.read()))
    if "ifc" not in config or "materiales" not in config or "hidratacion" not in config:
        raise ValueError("La configuración debe incluir ifc, materiales e hidratacion")
    return config


def aplicar_overrides(vertidas: list[VertidaIFC], config: dict) -> list[VertidaIFC]:
    sobrescritas = config.get("sobrescribir_vertidas", {})
    resultado = []
    for vertida in vertidas:
        cambio = sobrescritas.get(str(vertida.vertida_id), {}) or {}
        inicio = cambio.get("inicio_h")
        temperatura = cambio.get("temperatura_colocacion_c")
        resultado.append(replace(
            vertida,
            inicio_h=vertida.inicio_h if inicio is None else float(inicio),
            temperatura_colocacion_c=(vertida.temperatura_colocacion_c
                                      if temperatura is None else float(temperatura)),
        ))
    return resultado


def crear_fronteras(config: dict):
    borde = config.get("condicion_superficie", {"tipo": "aislado"})
    tipo = borde.get("tipo", "aislado").lower()
    if tipo == "aislado":
        condicion = aislado()
    elif tipo == "dirichlet":
        condicion = dirichlet(float(borde["temperatura_c"]))
    elif tipo == "neumann":
        condicion = neumann(float(borde.get("flujo_entrante_w_m2", 0.0)))
    elif tipo == "robin":
        condicion = robin(float(borde["coeficiente_h_w_m2_k"]),
                          float(borde["temperatura_aire_c"]))
    else:
        raise ValueError("tipo de condicion_superficie debe ser aislado, dirichlet, neumann o robin")
    return {cara: condicion for cara in ("x-", "x+", "y-", "y+", "z-", "z+")}


def campos_materiales(grilla, vertidas, materiales):
    """Construye k, rho y cp por celda a partir del tipo de cada vertida."""
    base = next(iter(materiales.values()))
    k = np.full(grilla.forma, float(base["conductividad_w_m_k"]))
    rho = np.full(grilla.forma, float(base["densidad_kg_m3"]))
    cp = np.full(grilla.forma, float(base["calor_especifico_j_kg_k"]))
    por_id = {v.vertida_id: v for v in vertidas}
    for identificador in np.unique(grilla.vertida_id):
        if identificador == 0:
            continue
        v = por_id[int(identificador)]
        if v.tipo_hormigon not in materiales:
            raise ValueError(f"Faltan propiedades térmicas para tipo {v.tipo_hormigon!r}")
        material = materiales[v.tipo_hormigon]
        seleccion = grilla.vertida_id == identificador
        k[seleccion] = float(material["conductividad_w_m_k"])
        rho[seleccion] = float(material["densidad_kg_m3"])
        cp[seleccion] = float(material["calor_especifico_j_kg_k"])
    return k, rho, cp


def guardar_convergencia_volumen(vertidas, resoluciones, ruta_csv):
    filas = []
    for dx in resoluciones:
        grilla = voxelizar_vertidas(vertidas, float(dx))
        filas.append({
            "dx_m": float(dx),
            "nx": grilla.forma[0], "ny": grilla.forma[1], "nz": grilla.forma[2],
            "voxeles_solidos": int(grilla.mascara_solido.sum()),
            "volumen_ifc_m3": sum(grilla.volumen_ifc_por_vertida_m3.values()),
            "volumen_voxelizado_m3": grilla.volumen_voxelizado_m3,
            "error_relativo_total": grilla.error_volumen_relativo_total,
        })
        print(f"  dx={dx:g} m: grilla={grilla.forma}, "
              f"Vvox={grilla.volumen_voxelizado_m3:.3f} m³, "
              f"error={100 * grilla.error_volumen_relativo_total:.3f}%")
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=filas[0].keys())
        escritor.writeheader()
        escritor.writerows(filas)
    return filas


def dibujar_geometria(grilla, vertidas_mostradas, ruta):
    ids = grilla.vertida_id
    visibles = (ids > 0) if vertidas_mostradas is None else np.isin(ids, vertidas_mostradas)
    colores = np.zeros(ids.shape + (4,), dtype=float)
    for identificador in np.unique(ids[visibles]):
        color = matplotlib.colors.to_rgba(COLORES_VERTIDA.get(int(identificador), "#777777"))
        colores[ids == identificador] = color
    fig = plt.figure(figsize=(9.2, 7.0))
    eje = fig.add_subplot(111, projection="3d")
    eje.voxels(visibles, facecolors=colores, edgecolor="#303030", linewidth=0.12)
    # Matplotlib dibuja los voxeles en coordenadas de índice (0, 1, ...).
    # Etiquetamos esas posiciones en metros para que la figura sea interpretable
    # junto con las dimensiones del IFC.
    for indice, eje_coord in enumerate((eje.xaxis, eje.yaxis, eje.zaxis)):
        cantidad = grilla.forma[indice]
        marcas = np.linspace(0, cantidad, min(6, cantidad + 1), dtype=int)
        valores_m = grilla.origen_m[indice] + marcas * grilla.paso_m[indice]
        eje_coord.set_ticks(marcas)
        eje_coord.set_ticklabels([f"{valor:.1f}" for valor in valores_m])
    eje.set_xlabel("x [m]")
    eje.set_ylabel("y [m]")
    eje.set_zlabel("z [m]")
    eje.set_title("Fundación del IFC voxelizada · color por vertida")
    eje.set_box_aspect((grilla.forma[0] * grilla.paso_m[0],
                        grilla.forma[1] * grilla.paso_m[1],
                        grilla.forma[2] * grilla.paso_m[2]))
    leyenda = [plt.Line2D([0], [0], marker="s", linestyle="", color=color,
                          label=f"Vertida {identificador}")
               for identificador, color in COLORES_VERTIDA.items()
               if np.any(ids[visibles] == identificador)]
    if leyenda:
        eje.legend(handles=leyenda, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.tight_layout()
    fig.savefig(ruta, dpi=170)
    return fig


def _clave_triangulo(vertices: np.ndarray, triangulo) -> tuple:
    """Clave geométrica independiente del orden de los vértices."""
    return tuple(sorted(
        tuple(np.round(vertices[int(indice)], 6)) for indice in triangulo
    ))


def dibujar_geometria_ifc(vertidas: list[VertidaIFC], ruta,
                          vertidas_mostradas=None):
    """Dibuja las mallas trianguladas que ya entrega ``extraer_vertidas``.

    No vuelve a leer el IFC ni crea un extractor paralelo: reutiliza
    ``VertidaIFC.vertices_m`` y ``VertidaIFC.triangulos``. Las caras internas
    que coinciden exactamente en las juntas entre vertidas se omiten del
    dibujo para evitar z-fighting; no se mueve ningún vértice ni cambia el
    volumen.
    """
    seleccionadas = [v for v in vertidas
                     if vertidas_mostradas is None or v.vertida_id in vertidas_mostradas]
    if not seleccionadas:
        raise ValueError("No hay vertidas para dibujar")

    apariciones = {}
    for vertida in seleccionadas:
        for triangulo in vertida.triangulos:
            clave = _clave_triangulo(vertida.vertices_m, triangulo)
            apariciones.setdefault(clave, []).append(vertida.vertida_id)
    claves_internas = {
        clave for clave, ids in apariciones.items() if len(set(ids)) > 1
    }

    fig = plt.figure(figsize=(10.6, 7.6))
    eje = fig.add_subplot(111, projection="3d")
    eje.set_proj_type("ortho")
    caras_todas = []
    colores_caras = []
    triangulos_omitidos = 0
    for vertida in seleccionadas:
        color = matplotlib.colors.to_rgba(
            COLORES_VERTIDA.get(vertida.vertida_id, "#777777")
        )
        for triangulo in vertida.triangulos:
            clave = _clave_triangulo(vertida.vertices_m, triangulo)
            if clave in claves_internas:
                triangulos_omitidos += 1
                continue
            caras_todas.append(vertida.vertices_m[triangulo])
            colores_caras.append(color)
    # Una sola colección permite que Matplotlib ordene todas las caras por
    # profundidad. Con una colección por vertida, las vertidas inferiores
    # pueden ocultar a las superiores por el orden de dibujo entre artistas.
    coleccion = Poly3DCollection(
        caras_todas, facecolors=colores_caras, edgecolors="#263442",
        linewidths=0.16, alpha=1.0, zsort="average",
    )
    eje.add_collection3d(coleccion)

    minimo = np.min(np.stack([v.minimo_m for v in seleccionadas]), axis=0)
    maximo = np.max(np.stack([v.maximo_m for v in seleccionadas]), axis=0)
    margen = 0.035 * np.maximum(maximo - minimo, 1.0)
    eje.set_xlim(minimo[0] - margen[0], maximo[0] + margen[0])
    eje.set_ylim(minimo[1] - margen[1], maximo[1] + margen[1])
    eje.set_zlim(minimo[2] - margen[2], maximo[2] + margen[2])
    eje.set_box_aspect(tuple(maximo - minimo), zoom=1.04)
    for eje_coord in (eje.xaxis, eje.yaxis, eje.zaxis):
        eje_coord.set_major_locator(MaxNLocator(5))
        eje_coord.set_tick_params(labelsize=8, pad=0)
    eje.set_xlabel("x [m]", labelpad=4)
    eje.set_ylabel("y [m]", labelpad=4)
    eje.set_zlabel("z [m]", labelpad=4)
    eje.view_init(elev=28, azim=-55)
    eje.set_title(
        "Geometría BIM de la fundación · mallas trianguladas del IFC\n"
        "Vista ortográfica con proporciones reales",
        fontsize=11.5, pad=10,
    )

    leyenda = [
        Patch(
            facecolor=COLORES_VERTIDA.get(v.vertida_id, "#777777"),
            edgecolor="#263442", alpha=1.0,
            label=(f"V{v.vertida_id} · inicio {v.inicio_h:g} h\n"
                   f"{v.volumen_m3:.2f} m³ · {len(v.triangulos)} triángulos IFC"),
        )
        for v in seleccionadas
    ]
    eje.legend(
        handles=leyenda, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        fontsize=7.8, frameon=True, title="Vertidas (mallas IFC)",
        title_fontsize=8.5, labelspacing=0.9, borderpad=0.7,
    )
    volumen = sum(v.volumen_m3 for v in seleccionadas)
    total_triangulos = sum(len(v.triangulos) for v in seleccionadas)
    detalle = (
        f"Dibujo: {total_triangulos - triangulos_omitidos} de {total_triangulos} "
        f"triángulos IFC; se omiten las caras internas coincidentes en las juntas.\n"
        if triangulos_omitidos else
        f"Dibujo: {total_triangulos} triángulos IFC; sin caras internas coincidentes.\n"
    )
    eje.text2D(0.015, 0.025,
               f"{len(seleccionadas)} mallas cerradas · {volumen:.3f} m³\n{detalle}"
               "Sin mover vértices ni alterar la geometría ni el volumen.",
               transform=eje.transAxes, fontsize=7.6, linespacing=1.35,
               bbox={"boxstyle": "round,pad=0.4", "facecolor": "white",
                     "edgecolor": "#dce3eb", "alpha": 0.95})
    fig.tight_layout()
    fig.savefig(ruta, dpi=180, bbox_inches="tight", pad_inches=0.15)
    return fig


def _normal_centroide(vertices: np.ndarray, triangulo):
    """Normal unitaria y centroide de un triángulo de la malla."""
    a, b, c = vertices[list(triangulo)]
    normal = np.cross(b - a, c - a)
    norma = float(np.linalg.norm(normal))
    if norma == 0.0:
        return np.zeros(3), (a + b + c) / 3.0
    return normal / norma, (a + b + c) / 3.0


def _envolvente_convexa(puntos: np.ndarray) -> np.ndarray:
    """Envolvente convexa 2D (monotone chain), sin dependencias externas."""
    puntos = sorted({tuple(np.round(p, 6)) for p in np.asarray(puntos)})
    if len(puntos) <= 2:
        return np.array(puntos, dtype=float)

    def cruz(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    inferior = []
    for punto in puntos:
        while len(inferior) >= 2 and cruz(inferior[-2], inferior[-1], punto) <= 1e-12:
            inferior.pop()
        inferior.append(punto)
    superior = []
    for punto in reversed(puntos):
        while len(superior) >= 2 and cruz(superior[-2], superior[-1], punto) <= 1e-12:
            superior.pop()
        superior.append(punto)
    return np.array(inferior[:-1] + superior[:-1], dtype=float)


def _huella_horizontal(vertida: VertidaIFC, z: float) -> np.ndarray:
    """Contorno 2D de la cara horizontal (tapa o base) en ``z``."""
    puntos = []
    for triangulo in vertida.triangulos:
        normal, centroide = _normal_centroide(vertida.vertices_m, triangulo)
        if (abs(abs(float(normal[2])) - 1.0) <= 1e-6
                and abs(float(centroide[2]) - z) <= 1e-6):
            puntos.extend(vertida.vertices_m[triangulo][:, :2])
    if not puntos:
        raise ValueError(
            f"La vertida V{vertida.vertida_id} no tiene cara horizontal en z={z}"
        )
    return _envolvente_convexa(np.asarray(puntos, dtype=float))


def _contiene_huella(exterior: np.ndarray, interior: np.ndarray) -> bool:
    """Indica si toda la huella ``interior`` está dentro de ``exterior``."""
    if len(exterior) < 3 or len(interior) == 0:
        return False
    for punto in interior:
        for i in range(len(exterior)):
            a = exterior[i]
            b = exterior[(i + 1) % len(exterior)]
            dx = float(b[0] - a[0])
            dy = float(b[1] - a[1])
            longitud = float(np.hypot(dx, dy))
            if longitud == 0.0:
                continue
            cruz = dx * (punto[1] - a[1]) - dy * (punto[0] - a[0])
            # Distancia firmada a la arista con tolerancia de 10 µm
            # (el IFC aproxima los planos del octágono a ~3 µm).
            if cruz / longitud < -1e-5:
                return False
    return True


def triangulos_exteriores_ensamblados(vertidas: list[VertidaIFC]):
    """Triángulos visibles de la envolvente exterior del ensamblaje.

    Solo se omiten las caras horizontales internas (tapas y bases) que quedan
    completamente cubiertas por otra vertida. Las caras laterales se conservan
    todas. No se modifica ningún vértice, dimensión ni volumen.
    """
    ordenadas = sorted(vertidas, key=lambda v: (v.minimo_m[2], v.maximo_m[2]))
    huellas = {}
    for vertida in ordenadas:
        huellas[(vertida.vertida_id, "top")] = _huella_horizontal(
            vertida, float(vertida.maximo_m[2])
        )
        huellas[(vertida.vertida_id, "bottom")] = _huella_horizontal(
            vertida, float(vertida.minimo_m[2])
        )

    visibles = {vertida.vertida_id: set() for vertida in ordenadas}
    for vertida in ordenadas:
        for indice, triangulo in enumerate(vertida.triangulos):
            normal, _ = _normal_centroide(vertida.vertices_m, triangulo)
            horizontal = abs(abs(float(normal[2])) - 1.0) <= 1e-6
            if not horizontal:
                visibles[vertida.vertida_id].add(indice)
                continue
            omitir = False
            if normal[2] > 0:
                for otra in ordenadas:
                    if otra is vertida:
                        continue
                    if abs(float(otra.minimo_m[2]) - float(vertida.maximo_m[2])) > 1e-6:
                        continue
                    if _contiene_huella(
                        huellas[(otra.vertida_id, "bottom")],
                        huellas[(vertida.vertida_id, "top")],
                    ):
                        omitir = True
                        break
            else:
                for otra in ordenadas:
                    if otra is vertida:
                        continue
                    if abs(float(otra.maximo_m[2]) - float(vertida.minimo_m[2])) > 1e-6:
                        continue
                    if _contiene_huella(
                        huellas[(otra.vertida_id, "top")],
                        huellas[(vertida.vertida_id, "bottom")],
                    ):
                        omitir = True
                        break
            if not omitir:
                visibles[vertida.vertida_id].add(indice)
    return visibles


def _aristas_caracteristicas(vertida: VertidaIFC, umbral_grados: float = 12.0):
    """Aristas donde cambia la orientación de las caras de la malla."""
    aristas = {}
    for triangulo in vertida.triangulos:
        normal, _ = _normal_centroide(vertida.vertices_m, triangulo)
        ids = [int(indice) for indice in triangulo]
        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            aristas.setdefault(tuple(sorted((a, b))), []).append(normal)
    coseno_umbral = float(np.cos(np.radians(umbral_grados)))
    resultado = []
    for (a, b), normales in aristas.items():
        if len(normales) == 1:
            resultado.append((a, b))
            continue
        destacar = False
        for i in range(len(normales)):
            for j in range(i + 1, len(normales)):
                if float(np.dot(normales[i], normales[j])) < coseno_umbral:
                    destacar = True
                    break
            if destacar:
                break
        if destacar:
            resultado.append((a, b))
    return resultado


def _caras_solidas(vertida: VertidaIFC, *, offset_z=0.0, triangulos=None):
    """Caras y colores sombreados de una vertida, sin añadirlos a la figura."""
    vertices = vertida.vertices_m.copy()
    vertices[:, 2] += offset_z
    indices = (range(len(vertida.triangulos)) if triangulos is None
               else sorted(triangulos))
    base = np.asarray(matplotlib.colors.to_rgb(
        COLORES_VERTIDA_VISTAS.get(vertida.vertida_id, "#777777")
    ))
    luz = np.asarray([-0.35, -0.55, 0.76])
    luz = luz / np.linalg.norm(luz)
    caras = []
    colores = []
    for indice in indices:
        triangulo = vertida.triangulos[indice]
        normal, _ = _normal_centroide(vertida.vertices_m, triangulo)
        factor = 0.70 + 0.30 * max(0.0, float(np.dot(normal, luz)))
        caras.append(vertices[triangulo])
        colores.append(tuple(np.clip(base * factor, 0.0, 1.0)) + (1.0,))
    return caras, colores


def _agregar_aristas(eje, vertida: VertidaIFC, *, offset_z=0.0,
                     ancho_arista=0.7):
    """Dibuja las aristas de faceta de una vertida sobre el sólido."""
    vertices = vertida.vertices_m.copy()
    vertices[:, 2] += offset_z
    for a, b in _aristas_caracteristicas(vertida):
        eje.plot(
            [vertices[a, 0], vertices[b, 0]],
            [vertices[a, 1], vertices[b, 1]],
            [vertices[a, 2], vertices[b, 2]],
            color="#111827", linewidth=ancho_arista,
            solid_capstyle="round", zorder=30,
        )


def _dibujar_solidos(eje, especificaciones):
    """Dibuja varias vertidas en una sola colección ordenada por profundidad.

    Usar una única ``Poly3DCollection`` permite que Matplotlib ordene las caras
    por profundidad y evita que una tapa grande de una capa inferior tape las
    capas superiores.
    """
    caras = []
    colores = []
    for especificacion in especificaciones:
        caras_vertida, colores_vertida = _caras_solidas(
            especificacion["vertida"],
            offset_z=especificacion.get("offset_z", 0.0),
            triangulos=especificacion.get("triangulos"),
        )
        caras.extend(caras_vertida)
        colores.extend(colores_vertida)

    coleccion = Poly3DCollection(
        caras, facecolors=colores, edgecolors="none", linewidths=0.0,
        alpha=1.0, zsort="average",
    )
    eje.add_collection3d(coleccion)

    for especificacion in especificaciones:
        if especificacion.get("mostrar_aristas", True):
            _agregar_aristas(
                eje, especificacion["vertida"],
                offset_z=especificacion.get("offset_z", 0.0),
                ancho_arista=especificacion.get("ancho_arista", 0.7),
            )
    return coleccion


def _dibujar_vertida_solida(eje, vertida: VertidaIFC, *, offset_z=0.0,
                            triangulos=None, mostrar_aristas=True,
                            ancho_arista=0.7):
    """Dibuja una vertida como sólido sombreado con contornos de faceta."""
    return _dibujar_solidos(eje, [{
        "vertida": vertida,
        "offset_z": offset_z,
        "triangulos": triangulos,
        "mostrar_aristas": mostrar_aristas,
        "ancho_arista": ancho_arista,
    }])


def _descripcion_geometrica(vertida: VertidaIFC) -> str:
    """Descripción breve de la forma de la vertida a partir de su malla."""
    z_inf = float(vertida.minimo_m[2])
    z_sup = float(vertida.maximo_m[2])
    anillo_inf = vertida.vertices_m[np.isclose(vertida.vertices_m[:, 2], z_inf, atol=1e-6)]
    anillo_sup = vertida.vertices_m[np.isclose(vertida.vertices_m[:, 2], z_sup, atol=1e-6)]
    ancho_inf = float(np.ptp(anillo_inf[:, 0]))
    ancho_sup = float(np.ptp(anillo_sup[:, 0]))
    if len(anillo_inf) > 12 or len(anillo_sup) > 12:
        return f"cilindro Ø{ancho_inf:.1f} m"
    if np.isclose(ancho_inf, ancho_sup):
        return f"prisma octogonal {ancho_inf:.1f}×{ancho_inf:.1f} m"
    return f"tronco octogonal {ancho_inf:.1f}→{ancho_sup:.1f} m"


def _frente_vertida(vertida: VertidaIFC, z: float) -> float:
    """Coordenada y de la cara frontal de la vertida a la altura ``z``."""
    z_inf = float(vertida.minimo_m[2])
    z_sup = float(vertida.maximo_m[2])
    if np.isclose(z_inf, z_sup):
        t = 0.0
    else:
        t = float(np.clip((z - z_inf) / (z_sup - z_inf), 0.0, 1.0))
    inferior = vertida.vertices_m[np.isclose(vertida.vertices_m[:, 2], z_inf, atol=1e-6)]
    superior = vertida.vertices_m[np.isclose(vertida.vertices_m[:, 2], z_sup, atol=1e-6)]
    return (1.0 - t) * float(inferior[:, 1].min()) + t * float(superior[:, 1].min())


def _configurar_vista_ifc(eje, *, minimo, maximo, zoom, elev, azim):
    """Ajusta límites, proporciones y cámara ortográfica de una vista 3D."""
    margen = 0.06 * np.maximum(maximo - minimo, 1.0)
    eje.set_xlim(minimo[0] - margen[0], maximo[0] + margen[0])
    eje.set_ylim(minimo[1] - margen[1], maximo[1] + margen[1])
    eje.set_zlim(minimo[2] - margen[2], maximo[2] + margen[2])
    eje.set_box_aspect(tuple(np.maximum(maximo - minimo, 0.1)), zoom=zoom)
    for eje_coord in (eje.xaxis, eje.yaxis, eje.zaxis):
        eje_coord.set_major_locator(MaxNLocator(5))
        eje_coord.set_tick_params(labelsize=8, pad=0)
        eje_coord.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        eje_coord.pane.set_edgecolor("#cbd5e1")
    eje.set_xlabel("x [m]", labelpad=4)
    eje.set_ylabel("y [m]", labelpad=4)
    eje.set_zlabel("z [m]", labelpad=4)
    eje.view_init(elev=elev, azim=azim)
    eje.grid(True)


def _etiqueta_frente(vertida: VertidaIFC, z_mid: float, offset_z: float = 0.0):
    x_etiqueta = {1: -2.2, 2: 2.2, 3: -1.8, 4: 1.2, 5: 0.0}.get(
        vertida.vertida_id, 0.0
    )
    return (
        x_etiqueta,
        _frente_vertida(vertida, z_mid) - 0.25,
        z_mid + offset_z,
    )


def _leyenda_vertidas(vertidas: list[VertidaIFC]):
    return [
        Patch(
            facecolor=COLORES_VERTIDA_VISTAS.get(v.vertida_id, "#777777"),
            edgecolor="#111827",
            label=(f"V{v.vertida_id} · inicio {v.inicio_h:g} h\n"
                   f"{_descripcion_geometrica(v)}\n"
                   f"{v.volumen_m3:.2f} m³"),
        )
        for v in vertidas
    ]


def dibujar_geometria_ifc_ensamblada(vertidas: list[VertidaIFC], ruta):
    """Vista A: fundación ensamblada en las coordenadas IFC originales.

    Dibuja la envolvente exterior con colores sólidos y contornos de faceta.
    Las tapas y bases internas cubiertas entre vertidas se omiten solo del
    dibujo para evitar superficies coplanares en conflicto; la geometría IFC
    no se altera.
    """
    if not vertidas:
        raise ValueError("No hay vertidas para dibujar")
    ordenadas = sorted(vertidas, key=lambda v: (v.minimo_m[2], v.maximo_m[2]))
    visibles = triangulos_exteriores_ensamblados(ordenadas)

    fig = plt.figure(figsize=(11.8, 8.2))
    eje = fig.add_subplot(111, projection="3d")
    eje.set_proj_type("ortho")
    _dibujar_solidos(eje, [
        {
            "vertida": vertida,
            "triangulos": visibles[vertida.vertida_id],
            "mostrar_aristas": True,
            "ancho_arista": 0.7,
        }
        for vertida in ordenadas
    ])

    minimo = np.min(np.stack([v.minimo_m for v in ordenadas]), axis=0)
    maximo = np.max(np.stack([v.maximo_m for v in ordenadas]), axis=0)
    _configurar_vista_ifc(eje, minimo=minimo, maximo=maximo, zoom=1.12,
                          elev=32, azim=-55)

    for vertida in ordenadas:
        z_mid = 0.5 * (float(vertida.minimo_m[2]) + float(vertida.maximo_m[2]))
        x, y, z = _etiqueta_frente(vertida, z_mid)
        eje.text(
            x, y, z, f"V{vertida.vertida_id} · inicio {vertida.inicio_h:g} h",
            fontsize=8.6, fontweight="bold", color="#0f172a",
            ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white",
                  "edgecolor": COLORES_VERTIDA_VISTAS.get(
                      vertida.vertida_id, "#777777"),
                  "alpha": 0.92, "linewidth": 0.8},
        )

    eje.set_title(
        "Fundación ensamblada · geometría IFC original\n"
        "Vista isométrica ortográfica · superficies sólidas y contornos de capa",
        fontsize=12.5, pad=12,
    )
    eje.legend(
        handles=_leyenda_vertidas(ordenadas), loc="upper left",
        bbox_to_anchor=(1.01, 1.0), fontsize=8.0, frameon=True,
        title="Vertidas · hora de inicio", title_fontsize=9.0,
        labelspacing=0.85, borderpad=0.8,
    )
    total = sum(len(visibles[v.vertida_id]) for v in ordenadas)
    total_malla = sum(len(v.triangulos) for v in ordenadas)
    volumen = sum(v.volumen_m3 for v in ordenadas)
    fig.text(
        0.01, 0.012,
        f"5 vertidas · {total} de {total_malla} triángulos dibujados · "
        f"volumen IFC {volumen:.3f} m³\n"
        "Se omiten solo las caras horizontales internas cubiertas entre "
        "vertidas. Vértices, dimensiones y volúmenes sin cambios.",
        fontsize=7.8, linespacing=1.4,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white",
              "edgecolor": "#dce3eb", "alpha": 0.95},
    )
    fig.tight_layout()
    fig.savefig(ruta, dpi=180, bbox_inches="tight", pad_inches=0.18)
    return fig


def _offsets_explotados(vertidas: list[VertidaIFC], separacion_m: float = 0.85):
    ordenadas = sorted(vertidas, key=lambda v: (v.minimo_m[2], v.maximo_m[2]))
    return {
        vertida.vertida_id: indice * float(separacion_m)
        for indice, vertida in enumerate(ordenadas)
    }


def dibujar_geometria_ifc_explotada(vertidas: list[VertidaIFC], ruta,
                                    separacion_m: float = 0.85):
    """Vista B: las cinco vertidas separadas verticalmente solo para el dibujo.

    Cada vertida se dibuja completa, con su geometría y dimensiones originales.
    El desplazamiento es exclusivamente visual y no representa la geometría de
    simulación.
    """
    if not vertidas:
        raise ValueError("No hay vertidas para dibujar")
    if separacion_m <= 0:
        raise ValueError("La separación debe ser positiva")
    ordenadas = sorted(vertidas, key=lambda v: (v.minimo_m[2], v.maximo_m[2]))
    offsets = _offsets_explotados(ordenadas, separacion_m)

    fig = plt.figure(figsize=(11.8, 9.4))
    eje = fig.add_subplot(111, projection="3d")
    eje.set_proj_type("ortho")
    _dibujar_solidos(eje, [
        {
            "vertida": vertida,
            "offset_z": offsets[vertida.vertida_id],
            "mostrar_aristas": True,
            "ancho_arista": 0.75,
        }
        for vertida in ordenadas
    ])

    for vertida in ordenadas:
        z_mid = 0.5 * (float(vertida.minimo_m[2]) + float(vertida.maximo_m[2]))
        x, y, z = _etiqueta_frente(
            vertida, z_mid, offsets[vertida.vertida_id]
        )
        eje.text(
            x, y, z, f"V{vertida.vertida_id} · inicio {vertida.inicio_h:g} h",
            fontsize=8.6, fontweight="bold", color="#0f172a",
            ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white",
                  "edgecolor": COLORES_VERTIDA_VISTAS.get(
                      vertida.vertida_id, "#777777"),
                  "alpha": 0.92, "linewidth": 0.8},
        )

    minimo = np.min(np.stack([v.minimo_m for v in ordenadas]), axis=0)
    maximo = np.max(np.stack([v.maximo_m for v in ordenadas]), axis=0).astype(float)
    maximo[2] += offsets[ordenadas[-1].vertida_id]
    _configurar_vista_ifc(eje, minimo=minimo, maximo=maximo, zoom=1.05,
                          elev=30, azim=-55)
    eje.set_title(
        "Fundación explotada · visualización explicativa\n"
        f"Las cinco vertidas separadas {separacion_m:g} m solo para el dibujo",
        fontsize=12.5, pad=12,
    )
    eje.legend(
        handles=_leyenda_vertidas(ordenadas), loc="upper left",
        bbox_to_anchor=(1.01, 1.0), fontsize=8.0, frameon=True,
        title="Vertidas · hora de inicio", title_fontsize=9.0,
        labelspacing=0.85, borderpad=0.8,
    )
    fig.text(
        0.01, 0.012,
        "Visualización explicativa, no geometría de simulación.\n"
        f"Los desplazamientos verticales ({separacion_m:g} m entre capas) son "
        "solo para el dibujo. Cada vertida conserva su geometría, dimensiones, "
        "volumen y malla IFC.",
        fontsize=7.9, linespacing=1.4,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#fff7ed",
              "edgecolor": "#f59e0b", "alpha": 0.97, "linewidth": 1.0},
    )
    fig.tight_layout()
    fig.savefig(ruta, dpi=180, bbox_inches="tight", pad_inches=0.18)
    return fig


def dibujar_temperatura(grilla, resultado: ResultadoEtapas, ruta):
    campo = resultado.temperaturas_c[-1]
    activa = resultado.mascaras_activas[-1]
    coordenadas = grilla.coordenadas_centros_m()
    ix, iy, iz = np.where(activa)
    temperaturas = campo[activa]
    fig = plt.figure(figsize=(9.2, 7.0))
    eje = fig.add_subplot(111, projection="3d")
    puntos = eje.scatter(coordenadas[0][ix], coordenadas[1][iy], coordenadas[2][iz],
                         c=temperaturas, cmap="inferno", s=9, marker="s")
    fig.colorbar(puntos, ax=eje, label="Temperatura [°C]", shrink=0.72)
    eje.set_xlabel("x [m]")
    eje.set_ylabel("y [m]")
    eje.set_zlabel("z [m]")
    eje.set_title(f"Temperatura de voxeles activos · t={resultado.tiempos_s[-1]/3600:.1f} h")
    eje.set_box_aspect((grilla.forma[0] * grilla.paso_m[0],
                        grilla.forma[1] * grilla.paso_m[1],
                        grilla.forma[2] * grilla.paso_m[2]))
    fig.tight_layout()
    fig.savefig(ruta, dpi=170)
    return fig


def dibujar_cortes_termicos(grilla, resultado: ResultadoEtapas,
                            ijk: tuple[int, int, int], tiempos_s, ruta):
    """Cortes ortogonales interiores en los instantes solicitados.

    Cada fila corresponde a un tiempo y cada columna a un plano x, y o z que
    pasa por ``ijk``. Se usan exclusivamente los campos guardados por el solver.
    """
    tiempos = list(dict.fromkeys(float(t) for t in tiempos_s))
    if not tiempos:
        raise ValueError("Se requiere al menos un tiempo para los cortes")
    cortes_por_tiempo = [
        (tiempo, cortes_por_punto(grilla, resultado, ijk, tiempo))
        for tiempo in tiempos
    ]
    valores = np.concatenate([
        corte.temperaturas_c[np.isfinite(corte.temperaturas_c)].ravel()
        for _, cortes in cortes_por_tiempo for corte in cortes
        if np.any(np.isfinite(corte.temperaturas_c))
    ])
    minimo = float(np.min(valores))
    maximo = float(np.max(valores))
    if np.isclose(minimo, maximo):
        maximo = minimo + 1.0
    normalizacion = matplotlib.colors.Normalize(vmin=minimo, vmax=maximo)
    mapa = matplotlib.colormaps["inferno"].with_extremes(bad="#d5dbe1")
    nombres = {"x": ("y [m]", "z [m]"), "y": ("x [m]", "z [m]"),
               "z": ("x [m]", "y [m]")}
    titulos = {
        "x": "Corte vertical YZ",
        "y": "Corte vertical XZ",
        "z": "Corte horizontal XY",
    }
    met = metricas_tmax(resultado)

    fig, ejes = plt.subplots(
        len(cortes_por_tiempo), 3, squeeze=False,
        figsize=(16.2, 5.2 * len(cortes_por_tiempo)),
    )
    fig.subplots_adjust(left=0.06, right=0.895, top=0.86, bottom=0.13,
                        hspace=0.55, wspace=0.16)
    for fila, (tiempo, cortes) in enumerate(cortes_por_tiempo):
        hora_ref = cortes[0].tiempo_h
        if np.isclose(hora_ref, met.t_pico_h, atol=1e-6):
            etiqueta_fila = (
                f"t = {hora_ref:g} h  ·  máximo térmico global "
                f"(Tmax = {met.tmax_pico_c:.2f} °C)"
            )
            color_fila = "#b2182b"
        elif np.isclose(hora_ref, resultado.tiempos_s[-1] / 3600.0, atol=1e-6):
            etiqueta_fila = (
                f"t = {hora_ref:g} h  ·  final de la simulación "
                f"(Tmax = {met.tmax_final_c:.2f} °C)"
            )
            color_fila = "#2166ac"
        else:
            etiqueta_fila = f"t = {hora_ref:g} h"
            color_fila = "#334155"
        posicion = ejes[fila, 0].get_position()
        fig.text(0.06, posicion.y1 + 0.045, etiqueta_fila,
                 ha="left", va="bottom", fontsize=11, fontweight="bold",
                 color=color_fila)

        for columna, corte in enumerate(cortes):
            eje = ejes[fila, columna]
            otros = [d for d in range(3) if d != "xyz".index(corte.eje)]
            extension = (
                grilla.origen_m[otros[0]],
                grilla.origen_m[otros[0]] + grilla.forma[otros[0]] * grilla.paso_m[otros[0]],
                grilla.origen_m[otros[1]],
                grilla.origen_m[otros[1]] + grilla.forma[otros[1]] * grilla.paso_m[otros[1]],
            )
            eje.imshow(
                np.ma.masked_invalid(corte.temperaturas_c.T),
                origin="lower", extent=extension, aspect="equal",
                cmap=mapa, norm=normalizacion, interpolation="nearest",
            )
            x_punto = grilla.origen_m[otros[0]] + (ijk[otros[0]] + 0.5) * grilla.paso_m[otros[0]]
            y_punto = grilla.origen_m[otros[1]] + (ijk[otros[1]] + 0.5) * grilla.paso_m[otros[1]]
            eje.plot([x_punto], [y_punto], marker="+", color="white",
                     markeredgewidth=1.8, markersize=13, linestyle="")
            eje.plot([x_punto], [y_punto], marker="+", color="#111111",
                     markeredgewidth=0.7, markersize=13, linestyle="")
            eje.set_xlabel(nombres[corte.eje][0], fontsize=9.5)
            eje.set_ylabel(nombres[corte.eje][1], fontsize=9.5)
            eje.tick_params(labelsize=8)
            # Poda los ticks extremos para que las etiquetas de las esquinas
            # de x e y no se superpongan en los paneles horizontales.
            eje.xaxis.set_major_locator(MaxNLocator(5, prune="both"))
            eje.yaxis.set_major_locator(MaxNLocator(5, prune="both"))
            eje.set_title(
                f"{titulos[corte.eje]}  ·  {corte.eje} = "
                f"{corte.coordenada_m:.2f} m",
                fontsize=10.5, pad=6,
            )

    mapeador = matplotlib.cm.ScalarMappable(norm=normalizacion, cmap=mapa)
    mapeador.set_array([])
    barra_eje = fig.add_axes([0.915, 0.16, 0.013, 0.66])
    barra = fig.colorbar(mapeador, cax=barra_eje)
    barra.set_label("Temperatura del hormigón [°C]", fontsize=9.5)
    barra.ax.tick_params(labelsize=8)

    fig.suptitle("Cortes térmicos interiores de la fundación",
                 fontsize=15, fontweight="bold", y=0.985)
    fig.text(0.5, 0.945,
             "Planos ortogonales que pasan por P1, el voxel de máximo global",
             ha="center", fontsize=10, color="#475569")

    xyz = tuple(
        grilla.origen_m[d] + (ijk[d] + 0.5) * grilla.paso_m[d]
        for d in range(3)
    )
    explicacion = (
        f"Cruz (+): P1, voxel de máximo global (ijk = {ijk[0]},{ijk[1]},{ijk[2]}; "
        f"x = {xyz[0]:.2f}, y = {xyz[1]:.2f}, z = {xyz[2]:.2f} m).\n"
        "Gris: celdas exteriores al hormigón (sin temperatura).  "
        "Color: temperatura del hormigón en el plano [°C]."
    )
    fig.text(0.06, 0.025, explicacion, ha="left", va="bottom",
             fontsize=9.2, color="#334155", linespacing=1.45)
    fig.savefig(ruta, dpi=170, bbox_inches="tight")
    return fig


def dibujar_evolucion_temperatura(grilla, resultado: ResultadoEtapas,
                                 tiempos_muestreo_s, ruta):
    """Guarda varios campos 3D en una sola figura con escala térmica común."""
    indices = [int(np.argmin(np.abs(resultado.tiempos_s - tiempo)))
               for tiempo in tiempos_muestreo_s]
    campos = [resultado.temperaturas_c[indice][resultado.mascaras_activas[indice]]
              for indice in indices]
    valores = np.concatenate([campo for campo in campos if campo.size])
    minimo = float(np.min(valores))
    maximo = float(np.max(valores))
    if np.isclose(minimo, maximo):
        maximo = minimo + 1.0
    normalizacion = matplotlib.colors.Normalize(vmin=minimo, vmax=maximo)
    colores_mapa = matplotlib.colormaps["inferno"]
    coordenadas = grilla.coordenadas_centros_m()

    cantidad = len(indices)
    columnas = 3 if cantidad <= 6 else 4
    filas = int(np.ceil(cantidad / columnas))
    fig = plt.figure(figsize=(18.0 if columnas == 3 else 20.0, 5.0 * filas))
    ejes = []
    aspecto = tuple(grilla.forma[i] * grilla.paso_m[i] for i in range(3))
    for panel, indice in enumerate(indices, start=1):
        eje = fig.add_subplot(filas, columnas, panel, projection="3d")
        ejes.append(eje)
        mascara = resultado.mascaras_activas[indice]
        ix, iy, iz = np.where(mascara)
        temperatura = resultado.temperaturas_c[indice][mascara]
        eje.scatter(coordenadas[0][ix], coordenadas[1][iy], coordenadas[2][iz],
                    c=temperatura, cmap=colores_mapa, norm=normalizacion,
                    s=3, marker="s", linewidths=0)
        eje.set_title(f"t = {resultado.tiempos_s[indice] / 3600.0:g} h", fontsize=10, pad=1)
        for dimension, setter_limite, setter_marcas in (
                (0, eje.set_xlim, eje.set_xticks),
                (1, eje.set_ylim, eje.set_yticks),
                (2, eje.set_zlim, eje.set_zticks)):
            inicio = grilla.origen_m[dimension]
            fin = inicio + grilla.forma[dimension] * grilla.paso_m[dimension]
            setter_limite(inicio, fin)
            setter_marcas(np.linspace(inicio, fin, 3))
        eje.set_xlabel("x [m]", fontsize=8, labelpad=-2)
        eje.set_ylabel("y [m]", fontsize=8, labelpad=-2)
        eje.set_zlabel("z [m]", fontsize=8, labelpad=-2)
        # En paneles 3D pequeños, los números de x e y chocan en las esquinas;
        # se conservan las etiquetas de eje y las marcas de altura z.
        eje.set_xticklabels([])
        eje.set_yticklabels([])
        eje.tick_params(labelsize=7, pad=0)
        eje.set_box_aspect(aspecto)
        eje.view_init(elev=25, azim=-58)

    mapeador = matplotlib.cm.ScalarMappable(norm=normalizacion, cmap=colores_mapa)
    mapeador.set_array([])
    barra_eje = fig.add_axes([0.92, 0.22, 0.016, 0.56])
    barra = fig.colorbar(mapeador, cax=barra_eje)
    barra.set_label("Temperatura [°C]", fontsize=10)
    barra.ax.tick_params(labelsize=8)
    fig.suptitle("Evolución espacial de la temperatura en la fundación", y=0.97)
    fig.subplots_adjust(left=0.015, right=0.90, bottom=0.04, top=0.91,
                        wspace=0.03, hspace=0.02)
    fig.savefig(ruta, dpi=170)
    return fig


def guardar_historial(resultado: ResultadoEtapas, ruta_csv):
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(["tiempo_h", "vertidas_activadas_en_instante", "Tmax_activa_C", "Tmedia_activa_C"])
        for t, ids, maxima, media in zip(
                resultado.tiempos_s, resultado.ids_activados,
                resultado.temperaturas_maximas_c, resultado.temperaturas_promedio_c):
            escritor.writerow([f"{t / 3600:.6f}", ",".join(map(str, ids)),
                               f"{maxima:.8f}", f"{media:.8f}"])


def dibujar_historial(resultado: ResultadoEtapas, ruta):
    """Curva T_max(t) con pico y activaciones de vertidas explícitos."""
    horas = resultado.tiempos_s / 3600.0
    maximas = np.asarray(resultado.temperaturas_maximas_c, dtype=float)
    met = metricas_tmax(resultado)
    fig, eje = plt.subplots(figsize=(9.2, 5.2))
    eje.plot(horas, maximas, color="#e87522", marker="o", ms=3.5,
             linewidth=2.2, label="Máxima en hormigón activo")
    eje.plot(horas, resultado.temperaturas_promedio_c, color="#197f8e",
             marker="s", ms=3.0, linewidth=1.6,
             label="Media en hormigón activo")
    for hora, ids in eventos_activacion(resultado):
        eje.axvline(hora, color="#94a3b8", linestyle=":", linewidth=0.9, zorder=0)
        eje.text(hora, float(np.nanmax(maximas)), f"+V{','.join(map(str, ids))}",
                 rotation=90, va="top", ha="right", fontsize=7, color="#64748b")
    eje.scatter([met.t_pico_h], [met.tmax_pico_c], s=75, color="#c83d4d",
                edgecolor="white", linewidth=1.1, zorder=5,
                label=f"Pico {met.tmax_pico_c:.3f} °C @ {met.t_pico_h:.1f} h")
    eje.annotate(
        f"Tmax = {met.tmax_pico_c:.3f} °C\nt = {met.t_pico_h:.1f} h",
        xy=(met.t_pico_h, met.tmax_pico_c),
        xytext=(0.55, 0.12), textcoords="axes fraction",
        fontsize=9, color="#7f1d1d",
        arrowprops={"arrowstyle": "->", "color": "#c83d4d", "lw": 1.2},
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#fff7f5",
              "edgecolor": "#f1b7ad"},
    )
    eje.set_xlabel("Tiempo desde primera vertida [h]")
    eje.set_ylabel("Temperatura [°C]")
    eje.set_title("T_max(t) y temperatura media del hormigón activo")
    eje.grid(True, alpha=0.3)
    eje.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(ruta, dpi=170)
    return fig


def dibujar_puntos_representativos(grilla, resultado: ResultadoEtapas,
                                   puntos: tuple[PuntoAnalisis, ...], ruta):
    """Curvas T(t) de los tres puntos elegidos por ``seleccionar_puntos...``.

    Los máximos se resumen en una tabla lateral para no superponer anotaciones
    sobre las curvas. El panel derecho también explica el criterio de selección.
    """
    if not puntos:
        raise ValueError("Se requiere al menos un punto representativo")
    colores = ["#d94801", "#2171b5", "#238b45", "#6a51a3", "#cb181d"]
    nombres = {
        "pico_global": "Punto más caliente",
        "superficie_expuesta": "Superficie expuesta",
        "nucleo_ultima_vertida": "Núcleo de la última vertida",
    }
    criterios = {
        "pico_global": ("Voxel con la mayor temperatura\n"
                        "registrada en toda la simulación\n(V4, interior)."),
        "superficie_expuesta": ("Voxel de una cara exterior con\n"
                                "la mayor temperatura de su\n"
                                "historia (V3)."),
        "nucleo_ultima_vertida": ("Voxel interior de la última\n"
                                  "vertida (V5) más cercano a su\n"
                                  "centroide."),
    }

    registros = []
    for posicion, punto in enumerate(puntos):
        horas, temperaturas = serie_punto(resultado, punto.ijk)
        validos = np.isfinite(temperaturas)
        if not np.any(validos):
            continue
        indices = np.flatnonzero(validos)
        indice_pico = int(indices[np.argmax(temperaturas[validos])])
        registros.append({
            "punto": punto,
            "horas": horas,
            "temperaturas": temperaturas,
            "color": colores[posicion % len(colores)],
            "nombre": nombres.get(
                punto.categoria, punto.categoria.replace("_", " ").capitalize()
            ),
            "criterio": criterios.get(
                punto.categoria, "Punto representativo de la grilla."
            ),
            "pico": float(temperaturas[indice_pico]),
            "hora_pico": float(horas[indice_pico]),
            "final": float(temperaturas[-1]),
        })
    if not registros:
        raise ValueError("Ningún punto representativo tiene datos activos")

    fig = plt.figure(figsize=(13.8, 7.6))
    grilla_ejes = fig.add_gridspec(
        1, 2, width_ratios=[2.35, 1.0], wspace=0.12,
        left=0.07, right=0.985, top=0.90, bottom=0.12,
    )
    eje = fig.add_subplot(grilla_ejes[0, 0])
    panel = fig.add_subplot(grilla_ejes[0, 1])
    panel.set_axis_off()

    for registro in registros:
        punto = registro["punto"]
        eje.plot(
            registro["horas"], registro["temperaturas"],
            color=registro["color"], linewidth=2.3,
            label=(f"{punto.etiqueta} · {registro['nombre']} "
                   f"(V{punto.vertida_id})"),
        )
        eje.scatter(
            [registro["hora_pico"]], [registro["pico"]],
            color=registro["color"], edgecolor="white", linewidth=1.0,
            s=70, zorder=6,
        )

    for hora, ids in eventos_activacion(resultado):
        eje.axvline(hora, color="#94a3b8", linestyle=":", linewidth=0.9, zorder=0)
        eje.text(hora, 0.015, f"+V{','.join(map(str, ids))}", rotation=90,
                 transform=eje.get_xaxis_transform(), va="bottom", ha="center",
                 fontsize=7.5, color="#64748b")

    valores = np.concatenate([
        registro["temperaturas"][np.isfinite(registro["temperaturas"])]
        for registro in registros
    ])
    hora_final = float(resultado.tiempos_s[-1] / 3600.0)
    eje.set_xlim(-3.0, hora_final * 1.01)
    eje.set_ylim(np.floor(valores.min() - 1.0), np.ceil(valores.max() + 1.0))
    eje.set_xticks(np.arange(0.0, hora_final + 1e-9, 24.0))
    eje.set_xlabel("Tiempo desde la primera vertida [h]")
    eje.set_ylabel("Temperatura [°C]")
    eje.grid(True, alpha=0.25)
    eje.tick_params(labelsize=9)
    for lado in ("top", "right"):
        eje.spines[lado].set_visible(False)

    panel.text(0.0, 1.0, "Criterio de selección", transform=panel.transAxes,
               ha="left", va="top", fontsize=11.5, fontweight="bold",
               color="#0f172a")
    y_texto = 0.91
    for registro in registros:
        punto = registro["punto"]
        panel.text(
            0.0, y_texto,
            f"● {punto.etiqueta} · {registro['nombre']} (V{punto.vertida_id})",
            transform=panel.transAxes, ha="left", va="top",
            fontsize=9.2, fontweight="bold", color=registro["color"],
        )
        panel.text(
            0.0, y_texto - 0.045, registro["criterio"],
            transform=panel.transAxes, ha="left", va="top",
            fontsize=8.2, color="#334155", linespacing=1.25,
        )
        y_texto -= 0.17

    panel.text(0.0, 0.37, "Máximos (● en la curva) y valor final",
               transform=panel.transAxes, ha="left", va="top",
               fontsize=9.2, fontweight="bold", color="#0f172a")
    filas = [
        [registro["punto"].etiqueta,
         f"{registro['pico']:.2f}",
         f"{registro['hora_pico']:.0f}",
         f"{registro['final']:.2f}"]
        for registro in registros
    ]
    tabla = panel.table(
        cellText=filas,
        colLabels=["Punto", "Pico [°C]", "t [h]", "Final [°C]"],
        cellLoc="center", colLoc="center",
        bbox=[0.0, 0.015, 1.0, 0.29],
    )
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(8.2)
    for (fila_tabla, _), celda in tabla.get_celld().items():
        celda.set_edgecolor("#cbd5e1")
        celda.set_linewidth(0.6)
        if fila_tabla == 0:
            celda.set_facecolor("#e2e8f0")
            celda.set_text_props(fontweight="bold", color="#0f172a")
        else:
            celda.set_facecolor("#ffffff")

    fig.suptitle("Temperaturas en tres puntos representativos",
                 fontsize=15, fontweight="bold", y=0.975)
    fig.text(0.07, 0.925,
             "Fundación del IFC · curvas del mismo campo térmico 3D del solver; "
             "líneas punteadas = colocación de cada vertida.",
             ha="left", fontsize=9, color="#475569")
    fig.savefig(ruta, dpi=170, bbox_inches="tight")
    return fig


def calcular_perfil_vertical_nucleo_superficie(grilla, resultado: ResultadoEtapas,
                                               ijk_referencia=None):
    """Diferencia vertical núcleo-superficie superior en una columna.

    Método de posprocesamiento, sin tocar el solver ni los datos térmicos. Se
    elige la columna (i, j) del voxel más caliente de toda la simulación. En
    cada instante se identifica la superficie superior como el voxel activo de
    mayor ``k`` y el núcleo como el voxel interior más caliente de esa misma
    columna, excluyendo los extremos superior e inferior. La métrica es
    ``T_núcleo - T_superficie`` y es la más comparable con la diferencia
    núcleo-superficie vertical reportada por Ju y Lei (2019).
    """
    if not resultado.temperaturas_c:
        raise ValueError("El resultado no contiene estados térmicos")
    if ijk_referencia is None:
        indice_pico = int(np.argmax(resultado.temperaturas_maximas_c))
        campo_pico = resultado.temperaturas_c[indice_pico]
        activa_pico = resultado.mascaras_activas[indice_pico]
        ijk_referencia = np.unravel_index(
            int(np.argmax(np.where(activa_pico, campo_pico, -np.inf))),
            campo_pico.shape,
        )
    i, j = int(ijk_referencia[0]), int(ijk_referencia[1])
    ejes = grilla.coordenadas_centros_m()
    filas = []
    for tiempo_s, campo, activa in zip(
            resultado.tiempos_s, resultado.temperaturas_c,
            resultado.mascaras_activas):
        ks = np.flatnonzero(activa[i, j, :])
        if ks.size < 2:
            continue
        k_superficie = int(ks[-1])
        k_inferior = int(ks[0])
        interiores = ks[(ks != k_superficie) & (ks != k_inferior)]
        if interiores.size == 0:
            continue
        k_nucleo = int(interiores[np.argmax(campo[i, j, interiores])])
        t_nucleo = float(campo[i, j, k_nucleo])
        t_superficie = float(campo[i, j, k_superficie])
        filas.append((
            float(tiempo_s) / 3600.0,
            float(ejes[0][i]), float(ejes[1][j]),
            float(ejes[2][k_nucleo]), t_nucleo,
            float(ejes[2][k_superficie]), t_superficie,
            t_nucleo - t_superficie,
        ))
    if not filas:
        raise ValueError("No hay columna vertical con núcleo y superficie")
    return filas, (i, j)


def guardar_perfil_vertical_e09(ruta_csv, filas) -> None:
    """Exporta la columna usada para la diferencia núcleo-superficie."""
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow([
            "tiempo_h", "x_m", "y_m", "z_nucleo_m", "t_nucleo_c",
            "z_superficie_m", "t_superficie_c", "delta_c",
        ])
        for fila in filas:
            escritor.writerow([f"{valor:.8f}" for valor in fila])


def guardar_metricas_e09(ruta_csv, *, grilla, resultado: ResultadoEtapas,
                         puntos: tuple[PuntoAnalisis, ...], hidratacion,
                         config: dict, dx: float, area_exterior_m2: float,
                         interfaces: dict,
                         perfil_vertical=None) -> None:
    """Resume en CSV las métricas que sustentan el reporte E09."""
    met = metricas_tmax(resultado)
    activaciones = "; ".join(
        f"{hora:g}h:V{','.join(map(str, ids))}"
        for hora, ids in eventos_activacion(resultado)
    )
    superficie = config.get("condicion_superficie", {})
    filas = [
        ("dx_m", f"{dx:.6f}", "m"),
        ("forma_grilla", "x".join(map(str, grilla.forma)), "voxeles"),
        ("voxeles_solidos", str(int(grilla.mascara_solido.sum())), "voxeles"),
        ("volumen_ifc_m3", f"{sum(grilla.volumen_ifc_por_vertida_m3.values()):.6f}", "m3"),
        ("volumen_voxelizado_m3", f"{grilla.volumen_voxelizado_m3:.6f}", "m3"),
        ("error_volumen_relativo", f"{grilla.error_volumen_relativo_total:.8%}", "-"),
        ("area_exterior_estimada_m2", f"{area_exterior_m2:.3f}", "m2"),
        ("interfaces_vertidas_m2", str(interfaces), "m2"),
        ("tmax_pico", f"{met.tmax_pico_c:.8f}", "C"),
        ("tmax_pico_h", f"{met.t_pico_h:.6f}", "h"),
        ("tmax_final", f"{met.tmax_final_c:.8f}", "C"),
        ("tmedia_final", f"{met.tmedia_final_c:.8f}", "C"),
        ("n_muestras_temporales", str(met.n_muestras), "-"),
        ("paso_maximo_estable", f"{resultado.ultimo_paso_maximo_estable_s:.6f}", "s"),
        ("residuo_energetico", f"{resultado.error_balance_j:.8e}", "J"),
        ("qmax", f"{hidratacion.qmax_w_m3:g}", "W/m3"),
        ("duracion_pulso", f"{hidratacion.duracion_s / 3600:g}", "h"),
        ("condicion_superficie", str(superficie.get("tipo", "aislado")), "-"),
        ("coef_convectivo_h", str(superficie.get("coeficiente_h_w_m2_k", "")), "W/m2K"),
        ("temperatura_aire", str(superficie.get("temperatura_aire_c", "")), "C"),
        ("activaciones", activaciones, "-"),
        ("puntos_representativos", ",".join(p.etiqueta for p in puntos), "-"),
    ]
    if perfil_vertical is not None:
        filas_perfil, (i_ref, j_ref) = perfil_vertical
        maximo = max(filas_perfil, key=lambda fila: fila[7])
        cercano_216 = min(filas_perfil, key=lambda fila: abs(fila[0] - 216.0))
        filas.extend([
            ("perfil_columna_ij", f"{i_ref},{j_ref}", "-"),
            ("dt_nucleo_superficie_max", f"{maximo[7]:.8f}", "C"),
            ("dt_nucleo_superficie_max_h", f"{maximo[0]:.6f}", "h"),
            ("dt_nucleo_superficie_216h", f"{cercano_216[7]:.8f}", "C"),
            ("dt_nucleo_superficie_216h_t_nucleo", f"{cercano_216[4]:.8f}", "C"),
            ("dt_nucleo_superficie_216h_t_superficie", f"{cercano_216[6]:.8f}", "C"),
        ])
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(["metrica", "valor", "unidad"])
        escritor.writerows(filas)


def guardar_puntos_e09(ruta_csv, resultado: ResultadoEtapas,
                       puntos: tuple[PuntoAnalisis, ...],
                       ruta_series=None) -> None:
    """Exporta el resumen de los puntos y, opcionalmente, sus curvas completas."""
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow([
            "etiqueta", "categoria", "ijk_x", "ijk_y", "ijk_z",
            "x_m", "y_m", "z_m", "vertida_id", "t_pico_h", "t_pico_c",
            "t_final_c", "n_muestras_activas",
        ])
        for punto in puntos:
            horas, temperaturas = serie_punto(resultado, punto.ijk)
            validos = np.isfinite(temperaturas)
            indice = int(np.flatnonzero(validos)[np.argmax(temperaturas[validos])])
            escritor.writerow([
                punto.etiqueta, punto.categoria, *punto.ijk,
                f"{punto.xyz_m[0]:.6f}", f"{punto.xyz_m[1]:.6f}",
                f"{punto.xyz_m[2]:.6f}", punto.vertida_id,
                f"{horas[indice]:.6f}", f"{temperaturas[indice]:.8f}",
                f"{temperaturas[-1]:.8f}", int(validos.sum()),
            ])

    if ruta_series is not None:
        with ruta_series.open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow(["tiempo_h", "etiqueta", "temperatura_c"])
            for punto in puntos:
                horas, temperaturas = serie_punto(resultado, punto.ijk)
                for hora, temperatura in zip(horas, temperaturas):
                    escritor.writerow([
                        f"{hora:.6f}", punto.etiqueta,
                        "" if not np.isfinite(temperatura) else f"{temperatura:.8f}",
                    ])



def calcular_control_termico(grilla, resultado: ResultadoEtapas,
                             tiempos_muestreo_s) -> dict[str, np.ndarray]:
    """Calcula Tmax activa y diferencia entre núcleo interior y superficie.

    La superficie se compone de voxeles con al menos una cara expuesta. El
    interior son los voxeles activos sin caras expuestas. La diferencia
    reportada es |Tmax(interior) - Tmin(superficie)|; queda como NaN mientras
    aún no exista un voxel interior.
    """
    indices = sorted({int(np.argmin(np.abs(resultado.tiempos_s - tiempo)))
                      for tiempo in tiempos_muestreo_s})
    tiempos = resultado.tiempos_s[indices]
    maximas = np.asarray(resultado.temperaturas_maximas_c[indices], dtype=float)
    maximas_interior = np.full(len(indices), np.nan, dtype=float)
    minimas_superficie = np.full(len(indices), np.nan, dtype=float)
    diferencias = np.full(len(indices), np.nan, dtype=float)

    for fila, indice in enumerate(indices):
        temperatura = resultado.temperaturas_c[indice]
        activa = resultado.mascaras_activas[indice]
        if not np.any(activa):
            continue
        caras = caras_superficie_exterior(grilla.vertida_id, activa)
        mascara_superficie = np.logical_or.reduce(tuple(caras.values()))
        mascara_interior = activa & ~mascara_superficie
        if not np.any(mascara_superficie) or not np.any(mascara_interior):
            continue
        t_interior_max = float(np.max(temperatura[mascara_interior]))
        t_superficie_min = float(np.min(temperatura[mascara_superficie]))
        maximas_interior[fila] = t_interior_max
        minimas_superficie[fila] = t_superficie_min
        diferencias[fila] = abs(t_interior_max - t_superficie_min)

    return {
        "tiempo_h": tiempos / 3600.0,
        "temperatura_maxima_c": maximas,
        "temperatura_maxima_interior_c": maximas_interior,
        "temperatura_minima_superficie_c": minimas_superficie,
        "diferencia_interior_superficie_c": diferencias,
    }


def guardar_control_termico(control: dict[str, np.ndarray], ruta_csv,
                            limite_temperatura_c: float,
                            limite_diferencia_c: float) -> None:
    """Guarda métricas térmicas y evaluación de límites en CSV."""
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow([
            "tiempo_h", "temperatura_maxima_c", "limite_temperatura_c",
            "excede_temperatura", "temperatura_maxima_interior_c",
            "temperatura_minima_superficie_c", "diferencia_interior_superficie_c",
            "limite_diferencia_c", "excede_diferencia",
        ])
        for indice, tiempo in enumerate(control["tiempo_h"]):
            tmax = float(control["temperatura_maxima_c"][indice])
            diferencia = float(control["diferencia_interior_superficie_c"][indice])
            diferencia_disponible = np.isfinite(diferencia)
            escritor.writerow([
                f"{tiempo:.6f}", f"{tmax:.8f}", f"{limite_temperatura_c:.6f}",
                tmax > limite_temperatura_c,
                _formato_control(control["temperatura_maxima_interior_c"][indice]),
                _formato_control(control["temperatura_minima_superficie_c"][indice]),
                _formato_control(diferencia), f"{limite_diferencia_c:.6f}",
                diferencia_disponible and diferencia > limite_diferencia_c,
            ])


def _formato_control(valor: float) -> str:
    return f"{float(valor):.8f}" if np.isfinite(valor) else ""


def caras_cubo_voxel(indice: tuple[int, int, int], origen_m, paso_m):
    """Devuelve los seis cuadriláteros de un voxel en coordenadas del modelo."""
    x0, y0, z0 = (origen_m[d] + indice[d] * paso_m[d] for d in range(3))
    x1, y1, z1 = x0 + paso_m[0], y0 + paso_m[1], z0 + paso_m[2]
    return [
        [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],
        [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
        [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
    ]


def dibujar_control_termico(control: dict[str, np.ndarray], ruta_png,
                            limite_temperatura_c: float,
                            limite_diferencia_c: float):
    """Grafica Tmax y diferencia interior-superficie con límites permitidos."""
    horas = control["tiempo_h"]
    tmax = control["temperatura_maxima_c"]
    diferencia = control["diferencia_interior_superficie_c"]
    fig, (eje_t, eje_delta) = plt.subplots(
        2, 1, figsize=(11.5, 7.6), sharex=True, facecolor="#f6f8fb",
        gridspec_kw={"hspace": 0.32},
    )
    for eje in (eje_t, eje_delta):
        eje.set_facecolor("white")
        eje.spines[["top", "right"]].set_visible(False)
        eje.grid(True, color="#dce3eb", linewidth=0.8, alpha=0.85)
        eje.set_axisbelow(True)

    eje_t.plot(horas, tmax, color="#e87522", marker="o", ms=4,
               linewidth=2.5, label="Máxima en todo el hormigón")
    eje_t.fill_between(horas, tmax, float(np.min(tmax)) - 2.0,
                       color="#e87522", alpha=0.09)
    eje_t.axhline(limite_temperatura_c, color="#c83d4d", linestyle=(0, (5, 3)),
                  linewidth=1.8, label=f"Límite: {limite_temperatura_c:g} °C")
    excede_t = tmax > limite_temperatura_c
    if np.any(excede_t):
        eje_t.scatter(horas[excede_t], tmax[excede_t], color="#c83d4d", s=38,
                      edgecolor="white", linewidth=0.7, zorder=4,
                      label="Sobre el límite")
    eje_t.set_ylabel("Temperatura [°C]")
    eje_t.set_title("Temperatura máxima del hormigón", loc="left", fontweight="bold")
    eje_t.set_ylim(float(np.min(tmax)) - 2.0,
                   max(limite_temperatura_c * 1.08, float(np.max(tmax)) + 2.0))
    eje_t.legend(loc="upper right", frameon=True, framealpha=0.95, edgecolor="#e3e8ef")

    eje_delta.plot(horas, diferencia, color="#197f8e", marker="o", ms=3.5,
                   linewidth=2.5, label="|T máx. interior − T mín. superficie|")
    diferencia_valida = diferencia[np.isfinite(diferencia)]
    if diferencia_valida.size:
        eje_delta.fill_between(horas, diferencia, 0.0, where=np.isfinite(diferencia),
                               color="#197f8e", alpha=0.09)
    eje_delta.axhline(limite_diferencia_c, color="#c83d4d", linestyle=(0, (5, 3)),
                      linewidth=1.8, label=f"Límite: {limite_diferencia_c:g} °C")
    excede_delta = np.isfinite(diferencia) & (diferencia > limite_diferencia_c)
    if np.any(excede_delta):
        eje_delta.scatter(horas[excede_delta], diferencia[excede_delta],
                          color="#c83d4d", s=38, edgecolor="white", linewidth=0.7,
                          zorder=4, label="Sobre el límite")
    eje_delta.set_xlabel("Tiempo desde la primera vertida [h]")
    eje_delta.set_ylabel("Diferencia térmica [°C]")
    eje_delta.set_title("Diferencia interior–superficie expuesta", loc="left",
                        fontweight="bold")
    eje_delta.set_ylim(0.0, max(limite_diferencia_c * 1.12,
                                float(np.max(diferencia_valida)) * 1.12
                                if diferencia_valida.size else 1.0))
    eje_delta.legend(loc="upper right", frameon=True, framealpha=0.95,
                     edgecolor="#e3e8ef")
    fig.suptitle("Control térmico de la fundación", fontsize=16,
                 fontweight="bold", color="#263442")
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.11, top=0.88, hspace=0.34)
    fig.savefig(ruta_png, dpi=190, facecolor=fig.get_facecolor())
    return fig


def crear_gif_control_termico_3d(grilla, resultado: ResultadoEtapas,
                                 control: dict[str, np.ndarray],
                                 limite_temperatura_c: float,
                                 limite_diferencia_c: float,
                                 ruta, fps: int = 2):
    """Anima curvas de límites y ubicación móvil de voxeles térmicos extremos."""
    indices = [int(np.argmin(np.abs(resultado.tiempos_s - tiempo * 3600.0)))
               for tiempo in control["tiempo_h"]]
    tiempos_h = np.asarray(control["tiempo_h"], dtype=float)
    aspecto = tuple(grilla.forma[d] * grilla.paso_m[d] for d in range(3))
    limites = [(grilla.origen_m[d],
                grilla.origen_m[d] + grilla.forma[d] * grilla.paso_m[d])
               for d in range(3)]
    nombres_cara = {"x-": 0, "x+": 1, "y-": 2, "y+": 3, "z-": 4, "z+": 5}
    mallas_por_estado = {}
    for indice in indices:
        activa = resultado.mascaras_activas[indice]
        ids_activos = tuple(sorted(int(i) for i in np.unique(grilla.vertida_id[activa])
                                   if i > 0))
        if ids_activos in mallas_por_estado:
            continue
        caras_expuestas = caras_superficie_exterior(grilla.vertida_id, activa)
        caras_malla = []
        for nombre, mascara_cara in caras_expuestas.items():
            for voxel in np.argwhere(mascara_cara):
                cubo = caras_cubo_voxel(tuple(int(i) for i in voxel),
                                        grilla.origen_m, grilla.paso_m)
                caras_malla.append(cubo[nombres_cara[nombre]])
        malla_estado = Poly3DCollection(
            caras_malla,
            facecolors=(0.78, 0.82, 0.86, 0.045),
            edgecolors=(0.40, 0.46, 0.52, 0.58),
            linewidths=0.38, zorder=1,
        )
        malla_estado.set_visible(False)
        mallas_por_estado[ids_activos] = malla_estado

    fig = plt.figure(figsize=(15.5, 8.0), facecolor="#f6f8fb")
    rejilla = fig.add_gridspec(2, 2, width_ratios=(0.92, 1.35),
                               left=0.055, right=0.96, top=0.95, bottom=0.09,
                               hspace=0.38, wspace=0.28)
    eje_tmax = fig.add_subplot(rejilla[0, 0])
    eje_delta = fig.add_subplot(rejilla[1, 0], sharex=eje_tmax)
    eje_3d = fig.add_subplot(rejilla[:, 1], projection="3d")
    eje_3d.computed_zorder = False

    for eje in (eje_tmax, eje_delta):
        eje.set_facecolor("white")
        eje.spines[["top", "right"]].set_visible(False)
        eje.grid(True, color="#dce3eb", linewidth=0.75, alpha=0.85)
        eje.set_axisbelow(True)
    eje_tmax.set_title("Temperatura máxima", loc="left", fontweight="bold")
    eje_tmax.set_ylabel("T [°C]")
    eje_tmax.set_xlim(float(tiempos_h[0]), max(1.0, float(tiempos_h[-1])))
    valores_tmax = control["temperatura_maxima_c"]
    eje_tmax.set_ylim(float(np.min(valores_tmax)) - 2.0,
                      max(limite_temperatura_c * 1.08,
                          float(np.max(valores_tmax)) + 2.0))
    eje_tmax.axhline(limite_temperatura_c, color="#c83d4d", linestyle="--",
                     linewidth=1.5, label=f"Límite {limite_temperatura_c:g} °C")
    linea_tmax, = eje_tmax.plot([], [], color="#e87522", linewidth=2.4,
                                marker="o", markersize=3.5, label="T máxima")
    punto_tmax, = eje_tmax.plot([], [], linestyle="", marker="o", markersize=8,
                                color="#e87522", markeredgecolor="white",
                                markeredgewidth=1.1, zorder=5)
    cursor_tmax = eje_tmax.axvline(tiempos_h[0], color="#64748b", linestyle=":")
    eje_tmax.legend(loc="upper right", fontsize=8, frameon=True)

    eje_delta.set_title("Diferencia interior–superficie", loc="left", fontweight="bold")
    eje_delta.set_xlabel("Tiempo [h]")
    eje_delta.set_ylabel("dT [°C]")
    eje_delta.set_xlim(float(tiempos_h[0]), max(1.0, float(tiempos_h[-1])))
    datos_delta = control["diferencia_interior_superficie_c"]
    validos_delta = datos_delta[np.isfinite(datos_delta)]
    eje_delta.set_ylim(0.0, max(limite_diferencia_c * 1.12,
                                float(np.max(validos_delta)) * 1.12
                                if validos_delta.size else 1.0))
    eje_delta.axhline(limite_diferencia_c, color="#c83d4d", linestyle="--",
                      linewidth=1.5, label=f"Límite {limite_diferencia_c:g} °C")
    linea_delta, = eje_delta.plot([], [], color="#197f8e", linewidth=2.4,
                                  marker="o", markersize=3.5,
                                  label="dT interior-superficie")
    punto_delta, = eje_delta.plot([], [], linestyle="", marker="o", markersize=8,
                                  color="#197f8e", markeredgecolor="white",
                                  markeredgewidth=1.1, zorder=5)
    cursor_delta = eje_delta.axvline(tiempos_h[0], color="#64748b", linestyle=":")
    eje_delta.legend(loc="upper right", fontsize=8, frameon=True)

    eje_3d.set_title("Vertidas activas · t = 0 h",
                     fontsize=12, pad=14, fontweight="bold")
    eje_3d.set_xlabel("x [m]")
    eje_3d.set_ylabel("y [m]")
    eje_3d.set_zlabel("z [m]")
    eje_3d.set_xlim(*limites[0])
    eje_3d.set_ylim(*limites[1])
    eje_3d.set_zlim(*limites[2])
    eje_3d.set_xticklabels([])
    eje_3d.set_yticklabels([])
    eje_3d.set_zticks(np.linspace(limites[2][0], limites[2][1], 3))
    eje_3d.tick_params(axis="z", labelsize=8, pad=0)
    eje_3d.set_box_aspect(aspecto)
    eje_3d.view_init(elev=25, azim=-58)
    for malla_estado in mallas_por_estado.values():
        eje_3d.add_collection3d(malla_estado)
    eje_3d.legend(handles=[
        plt.Line2D(
            [0], [0], marker="s", linestyle="", markersize=8,
            markerfacecolor="#4D7C8A", markeredgecolor="white",
            label="Voxel más frío activo",
        ),
        plt.Line2D(
            [0], [0], marker="s", linestyle="", markersize=8,
            markerfacecolor="#E76F51", markeredgecolor="white",
            label="Voxel más caliente activo",
        ),
    ], loc="upper left", fontsize=8, frameon=True)
    texto_extremos = eje_3d.text2D(
        0.02, 0.02, "", transform=eje_3d.transAxes, fontsize=8,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white",
              "edgecolor": "#dce3eb", "alpha": 0.92},
    )
    artistas_extremos = []

    def actualizar(frame):
        indice = indices[frame]
        tiempo_h = float(tiempos_h[frame])
        campo = resultado.temperaturas_c[indice]
        activa = resultado.mascaras_activas[indice]
        posiciones_activas = np.argwhere(activa)
        posiciones_calientes = posiciones_activas
        posiciones_frias = posiciones_activas
        caliente = tuple(posiciones_calientes[
            int(np.argmax(campo[activa]))])
        fria = tuple(posiciones_frias[
            int(np.argmin(campo[activa]))])
        ids_activos = tuple(sorted(int(i) for i in np.unique(grilla.vertida_id[activa])
                                   if i > 0))
        for estado, malla_estado in mallas_por_estado.items():
            malla_estado.set_visible(estado == ids_activos)

        for artista in artistas_extremos:
            artista.remove()
        artistas_extremos.clear()
        alpha_extremo = 0.70 if caliente == fria else 0.97
        for voxel, color in ((fria, "#4D7C8A"), (caliente, "#E76F51")):
            cubo_extremo = Poly3DCollection(
                caras_cubo_voxel(voxel, grilla.origen_m, grilla.paso_m),
                facecolors=color, edgecolors="white", linewidths=0.9,
                alpha=alpha_extremo, zorder=10,
            )
            eje_3d.add_collection3d(cubo_extremo)
            artistas_extremos.append(cubo_extremo)
        etiqueta_ids = ", ".join(map(str, ids_activos)) or "ninguna"
        eje_3d.set_title(f"Vertidas activas: {etiqueta_ids} · t = {tiempo_h:.0f} h",
                         fontsize=12, pad=14, fontweight="bold")
        texto_extremos.set_text(
            f"Caliente ijk=({caliente[0]+1},{caliente[1]+1},{caliente[2]+1})  "
            f"{campo[caliente]:.1f} °C\n"
            f"Frío ijk=({fria[0]+1},{fria[1]+1},{fria[2]+1})  {campo[fria]:.1f} °C"
        )

        linea_tmax.set_data(tiempos_h[:frame + 1], valores_tmax[:frame + 1])
        punto_tmax.set_data([tiempo_h], [valores_tmax[frame]])
        cursor_tmax.set_xdata([tiempo_h, tiempo_h])
        linea_delta.set_data(tiempos_h[:frame + 1], datos_delta[:frame + 1])
        delta_actual = datos_delta[frame]
        if np.isfinite(delta_actual):
            punto_delta.set_data([tiempo_h], [delta_actual])
        else:
            punto_delta.set_data([], [])
        cursor_delta.set_xdata([tiempo_h, tiempo_h])
        return linea_tmax, punto_tmax, linea_delta, punto_delta

    animacion = FuncAnimation(fig, actualizar, frames=len(indices), interval=500,
                              blit=False, repeat=True)
    animacion.save(ruta, writer=PillowWriter(fps=fps), dpi=90)
    plt.close(fig)
    return len(indices)


def imprimir_control_termico(control: dict[str, np.ndarray],
                             limite_temperatura_c: float,
                             limite_diferencia_c: float) -> str:
    """Imprime la evolución por intervalo y devuelve resumen de cumplimiento."""
    print("\nCONTROL TÉRMICO")
    print(" tiempo [h] | Tmax [°C] | diferencia T interior-superficie [°C] | estado")
    for indice, tiempo in enumerate(control["tiempo_h"]):
        tmax = float(control["temperatura_maxima_c"][indice])
        diferencia = float(control["diferencia_interior_superficie_c"][indice])
        estado = []
        if tmax > limite_temperatura_c:
            estado.append("EXCEDE Tmax")
        if np.isfinite(diferencia) and diferencia > limite_diferencia_c:
            estado.append("EXCEDE diferencia T")
        if not np.isfinite(diferencia):
            estado.append("dT no disponible: sin nucleo interior")
        print(f"{tiempo:10.2f} | {tmax:9.3f} | {_formato_control(diferencia):>27} | "
              f"{', '.join(estado) if estado else 'Dentro de límites'}")

    indice_tmax = int(np.argmax(control["temperatura_maxima_c"]))
    resumen = (
        f"Tmax pico={control['temperatura_maxima_c'][indice_tmax]:.3f} °C "
        f"a t={control['tiempo_h'][indice_tmax]:.2f} h "
        f"(límite {limite_temperatura_c:g} °C: "
        f"{'EXCEDIDO' if control['temperatura_maxima_c'][indice_tmax] > limite_temperatura_c else 'cumple'})."
    )
    validas = np.isfinite(control["diferencia_interior_superficie_c"])
    if np.any(validas):
        indices_validos = np.flatnonzero(validas)
        indice_delta = int(indices_validos[np.argmax(
            control["diferencia_interior_superficie_c"][validas])])
        delta_max = float(control["diferencia_interior_superficie_c"][indice_delta])
        resumen += (
            f" dT pico={delta_max:.3f} °C a t={control['tiempo_h'][indice_delta]:.2f} h "
            f"(límite {limite_diferencia_c:g} °C: "
            f"{'EXCEDIDO' if delta_max > limite_diferencia_c else 'cumple'})."
        )
    else:
        resumen += " dT interior-superficie no calculable: no hubo nucleo interior en las muestras."
    print("Resumen: " + resumen)
    return resumen


def crear_gif_evolucion_3d(grilla, resultado: ResultadoEtapas,
                           tiempos_muestreo_s, ruta, fps: int = 4):
    """Anima el campo 3D y las temperaturas central/promedio sincronizadas."""
    indices = [int(np.argmin(np.abs(resultado.tiempos_s - tiempo)))
               for tiempo in tiempos_muestreo_s]
    tiempos_h = resultado.tiempos_s[indices] / 3600.0
    coordenadas = grilla.coordenadas_centros_m()
    campos_activos = [resultado.temperaturas_c[i][resultado.mascaras_activas[i]]
                      for i in indices]
    temperaturas_todas = np.concatenate([campo for campo in campos_activos if campo.size])
    minimo = float(np.min(temperaturas_todas))
    maximo = float(np.max(temperaturas_todas))
    if np.isclose(minimo, maximo):
        maximo = minimo + 1.0
    normalizacion = matplotlib.colors.Normalize(vmin=minimo, vmax=maximo)
    mapa_colores = matplotlib.colormaps["inferno"]

    # Para una geometría irregular, "centro" es el voxel sólido más cercano
    # al centro de la caja envolvente del IFC.
    mascara_final = resultado.mascaras_activas[-1]
    ix, iy, iz = np.where(mascara_final)
    centro_caja = np.array([
        grilla.origen_m[d] + 0.5 * grilla.forma[d] * grilla.paso_m[d]
        for d in range(3)
    ])
    distancia2 = sum((coordenadas[d][(ix, iy, iz)[d]] - centro_caja[d]) ** 2
                     for d in range(3))
    posicion_central = (int(ix[np.argmin(distancia2)]),
                        int(iy[np.argmin(distancia2)]),
                        int(iz[np.argmin(distancia2)]))
    temperatura_centro = np.asarray([
        resultado.temperaturas_c[i][posicion_central]
        if resultado.mascaras_activas[i][posicion_central] else np.nan
        for i in indices
    ])
    temperatura_promedio = resultado.temperaturas_promedio_c[indices]
    valores_curvas = np.concatenate((temperatura_centro, temperatura_promedio))
    y_min = float(np.nanmin(valores_curvas))
    y_max = float(np.nanmax(valores_curvas))
    margen_y = max(1.0, 0.05 * (y_max - y_min))

    fig = plt.figure(figsize=(13.0, 6.2))
    rejilla_figura = fig.add_gridspec(1, 3, width_ratios=(1.05, 0.045, 1.0), wspace=0.48)
    eje_3d = fig.add_subplot(rejilla_figura[0, 0], projection="3d")
    eje_color = fig.add_subplot(rejilla_figura[0, 1])
    eje_tiempo = fig.add_subplot(rejilla_figura[0, 2])
    mapeador = matplotlib.cm.ScalarMappable(norm=normalizacion, cmap=mapa_colores)
    mapeador.set_array([])
    barra = fig.colorbar(mapeador, cax=eje_color)
    eje_color.set_title("°C", fontsize=9, pad=4)
    aspecto = tuple(grilla.forma[d] * grilla.paso_m[d] for d in range(3))
    limites = [(grilla.origen_m[d],
                grilla.origen_m[d] + grilla.forma[d] * grilla.paso_m[d])
               for d in range(3)]

    def actualizar(frame):
        indice = indices[frame]
        mascara = resultado.mascaras_activas[indice]
        campo = resultado.temperaturas_c[indice]
        x, y, z = np.where(mascara)
        eje_3d.clear()
        eje_3d.scatter(coordenadas[0][x], coordenadas[1][y], coordenadas[2][z],
                       c=campo[mascara], cmap=mapa_colores, norm=normalizacion,
                       s=7, marker="s", linewidths=0)
        eje_3d.set_title(f"Temperatura · t = {tiempos_h[frame]:.1f} h")
        eje_3d.set_xlabel("x [m]")
        eje_3d.set_ylabel("y [m]")
        eje_3d.set_zlabel("z [m]")
        eje_3d.set_xlim(*limites[0])
        eje_3d.set_ylim(*limites[1])
        eje_3d.set_zlim(*limites[2])
        eje_3d.set_xticklabels([])
        eje_3d.set_yticklabels([])
        eje_3d.set_zticks(np.linspace(limites[2][0], limites[2][1], 3))
        eje_3d.tick_params(axis="z", labelsize=7, pad=0)
        eje_3d.set_box_aspect(aspecto)
        eje_3d.view_init(elev=25, azim=-58)

        eje_tiempo.clear()
        eje_tiempo.plot(tiempos_h[:frame + 1], temperatura_centro[:frame + 1],
                        color="#76add5", linewidth=2.0, label="Centro")
        eje_tiempo.plot(tiempos_h[:frame + 1], temperatura_promedio[:frame + 1],
                        color="#f28e2b", linewidth=2.0, label="Promedio")
        eje_tiempo.scatter([tiempos_h[frame]], [temperatura_centro[frame]],
                           color="#76add5", s=20, zorder=3)
        eje_tiempo.scatter([tiempos_h[frame]], [temperatura_promedio[frame]],
                           color="#f28e2b", s=20, zorder=3)
        eje_tiempo.axvline(tiempos_h[frame], color="#555555", linestyle="--", linewidth=1.0)
        eje_tiempo.set_xlim(0.0, float(tiempos_h[-1]) if tiempos_h[-1] > 0 else 1.0)
        eje_tiempo.set_ylim(y_min - margen_y, y_max + margen_y)
        eje_tiempo.set_xlabel("Tiempo [h]")
        eje_tiempo.set_ylabel("Temperatura [°C]")
        eje_tiempo.set_title("Temperatura vs. tiempo")
        eje_tiempo.grid(True, alpha=0.3)
        eje_tiempo.legend(loc="upper left", fontsize=8)
        return []

    animacion = FuncAnimation(fig, actualizar, frames=len(indices),
                              interval=1000 / fps, blit=False, repeat=True)
    animacion.save(ruta, writer=PillowWriter(fps=fps), dpi=95)
    plt.close(fig)
    return len(indices)


def main() -> int:
    args = argumentos()
    modo_e09 = args.modo_e09
    salida = (args.salida or (SALIDA_E09 if modo_e09 else SALIDA)).expanduser().resolve()
    # En E09 las cinco figuras y los CSV se guardan con los resultados,
    # dentro de results/semana3_e09. Otros modos conservan su ruta histórica.
    evidencia_predeterminada = salida if modo_e09 else EVIDENCIA_E09
    evidencia = (args.evidencia or evidencia_predeterminada).expanduser().resolve()
    if modo_e09:
        archivo_salida = salida.parent / f"{salida.name}_archivo"
        movidos_salida = archivar_contenido(salida, archivo_salida)
        if movidos_salida:
            print(f"Salida previa archivada: {len(movidos_salida)} archivo(s) en "
                  f"{archivo_salida}")
        # Si se solicita una carpeta de evidencia distinta, archivarla solo
        # una vez. Evita la duplicación de respaldos cuando ambas coinciden.
        if evidencia != salida:
            archivo_evidencia = evidencia.parent / f"{evidencia.name}_archivo"
            movidos_evidencia = archivar_contenido(evidencia, archivo_evidencia)
            if movidos_evidencia:
                print(f"Evidencia previa archivada: {len(movidos_evidencia)} "
                      f"archivo(s) en {archivo_evidencia}")
    config = cargar_configuracion(args.config)
    ruta_ifc = (RAIZ / config["ifc"]).resolve()
    vertidas = aplicar_overrides(extraer_vertidas(ruta_ifc), config)
    print(f"IFC: {ruta_ifc}\nVertidas extraídas: {len(vertidas)}")
    for v in vertidas:
        print(f"  ID {v.vertida_id} {v.nombre}: inicio={v.inicio_h:g} h, "
              f"Tcol={v.temperatura_colocacion_c:g} °C, "
              f"dim={np.round(v.dimensiones_m, 3)} m, V={v.volumen_m3:.3f} m³")

    salida.mkdir(parents=True, exist_ok=True)
    evidencia.mkdir(parents=True, exist_ok=True)
    if modo_e09:
        print(f"Modo E09: salida exclusiva={salida}")
        print(f"Modo E09: evidencia exclusiva={evidencia}")
    resoluciones = config.get("resoluciones_error_volumen_m", [])
    if resoluciones:
        print("\nConvergencia geométrica de volumen:")
        guardar_convergencia_volumen(vertidas, resoluciones, salida / "error_volumen.csv")

    dx = float(args.dx if args.dx is not None else config["paso_voxel_m"])
    grilla = voxelizar_vertidas(vertidas, dx)
    interfaces = area_interfaces_entre_vertidas(grilla.vertida_id, grilla.paso_m)
    superficies = caras_superficie_exterior(grilla.vertida_id)
    area_exterior_m2 = sum(
        np.count_nonzero(superficies[cara])
        * float(np.prod([grilla.paso_m[i] for i in range(3) if i != "xyz".index(cara[0])]))
        for cara in superficies
    )
    print(f"\nVoxelización principal: dx={dx:g} m, forma={grilla.forma}, "
          f"voxeles de hormigón={int(grilla.mascara_solido.sum())}")
    print(f"Volumen IFC={sum(grilla.volumen_ifc_por_vertida_m3.values()):.3f} m³; "
          f"volumen voxelizado={grilla.volumen_voxelizado_m3:.3f} m³; "
          f"error={100 * grilla.error_volumen_relativo_total:.3f}%")
    print(f"Área externa estimada por caras de voxel: {area_exterior_m2:.2f} m²")
    print(f"Áreas de contacto entre vertidas [m²]: {interfaces}")

    materiales = config["materiales"]
    k, rho, cp = campos_materiales(grilla, vertidas, materiales)
    hcfg = config["hidratacion"]
    hidratacion = ModeloPulsoCosenoidal(
        float(hcfg["qmax_w_m3"]),
        float(hcfg["duracion_pulso_h"]) * 3600.0,
    )
    print(f"Fuente cosenoidal por vertida: Qmax={hidratacion.qmax_w_m3:g} W/m³, "
          f"duración={hidratacion.duracion_s / 3600:g} h, "
          f"pico a las {hidratacion.duracion_s / 7200:g} h de cada colocación.")

    final_s = float(args.horas if args.horas is not None else config["tiempo_final_h"]) * 3600.0
    control_cfg = config.get("control_termico", {})
    intervalo_control_h = float(control_cfg.get("intervalo_h", 6.0))
    if intervalo_control_h <= 0:
        raise ValueError("control_termico.intervalo_h debe ser positivo")
    intervalo_control_s = intervalo_control_h * 3600.0
    tiempos_control_s = np.arange(0.0, final_s + 1e-9, intervalo_control_s).tolist()
    if not tiempos_control_s:
        tiempos_control_s = [0.0]
    if not np.isclose(tiempos_control_s[-1], final_s):
        tiempos_control_s.append(final_s)

    # Guardar un perfil por día, calculado desde el tiempo final.
    # Así la configuración no necesita una lista manual de tiempos.
    intervalo_perfil_s = 24.0 * 3600.0
    tiempos_evolucion_s = np.arange(0.0, final_s + 1e-9, intervalo_perfil_s).tolist()
    if not tiempos_evolucion_s:
        tiempos_evolucion_s = [0.0]
    if not np.isclose(tiempos_evolucion_s[-1], final_s):
        tiempos_evolucion_s.append(final_s)
    intervalo_gif_h = float(config.get("intervalo_gif_h", 4.0))
    if intervalo_gif_h <= 0:
        raise ValueError("intervalo_gif_h debe ser positivo")
    intervalo_gif_s = intervalo_gif_h * 3600.0
    tiempos_gif_s = np.arange(0.0, final_s + 1e-9, intervalo_gif_s).tolist()
    if not tiempos_gif_s or not np.isclose(tiempos_gif_s[-1], final_s):
        tiempos_gif_s.append(final_s)
    tiempos_guardado_s = sorted(set(tiempos_evolucion_s + tiempos_gif_s + tiempos_control_s))
    resultado = simular_vertidas(
        grilla=grilla,
        vertidas=vertidas,
        conductividad_w_m_k=k,
        densidad_kg_m3=rho,
        calor_especifico_j_kg_k=cp,
        hidratacion=hidratacion,
        tiempo_final_s=final_s,
        fraccion_paso_estable=float(config.get("fraccion_paso_estable", 0.8)),
        paso_tiempo_max_s=float(config.get("paso_tiempo_max_s", 600.0)),
        tiempos_muestreo_s=tiempos_guardado_s,
        fronteras=crear_fronteras(config),
        progreso=lambda mensaje: print("  " + mensaje),
    )
    print(f"\nBalance energético: E_final={resultado.energia_final_j:.6e} J, "
          f"E_activación={resultado.energia_activacion_j:.6e} J, "
          f"E_fuente={resultado.energia_fuente_j:.6e} J, "
          f"E_frontera={resultado.energia_frontera_j:.6e} J, "
          f"residuo={resultado.error_balance_j:.6e} J")

    limite_temperatura_c = float(control_cfg.get("limite_temperatura_maxima_c", 70.0))
    limite_diferencia_c = float(control_cfg.get(
        "limite_diferencia_interior_superficie_c", 20.0))
    control = calcular_control_termico(grilla, resultado, tiempos_control_s)
    resumen_control = imprimir_control_termico(
        control, limite_temperatura_c, limite_diferencia_c,
    )
    guardar_control_termico(
        control, salida / "control_termico.csv",
        limite_temperatura_c, limite_diferencia_c,
    )

    (salida / "resumen.txt").write_text(
        "Pipeline IFC 3D por vertidas\n"
        f"IFC: {ruta_ifc}\n"
        f"dx: {dx:g} m; forma: {grilla.forma}; voxeles sólidos: {int(grilla.mascara_solido.sum())}\n"
        f"Volumen IFC: {sum(grilla.volumen_ifc_por_vertida_m3.values()):.6f} m3\n"
        f"Volumen voxelizado: {grilla.volumen_voxelizado_m3:.6f} m3\n"
        f"Error relativo total: {grilla.error_volumen_relativo_total:.6%}\n"
        f"Interfaces de vertidas [m2]: {interfaces}\n"
        f"Tmax final: {resultado.temperaturas_maximas_c[-1]:.6f} C\n"
        f"Tmedia final: {resultado.temperaturas_promedio_c[-1]:.6f} C\n"
        f"Residuo energético: {resultado.error_balance_j:.8e} J\n"
        f"Control térmico (cada {intervalo_control_h:g} h): {resumen_control}\n"
        f"Fuente cosenoidal: Qmax={hidratacion.qmax_w_m3:g} W/m3; "
        f"duración={hidratacion.duracion_s / 3600:g} h por vertida.\n"
        f"Condición superficial: {config.get('condicion_superficie', {}).get('tipo', 'aislado')}; "
        "misma condición en todas las caras exteriores.\n"
        "Simplificaciones: cinética dependiente de edad (no de temperatura), sin suelo/moldaje/clima variable, "
        "malla cartesiana staircase.\n",
        encoding="utf-8",
    )
    guardar_historial(resultado, salida / "historial_termico.csv")

    ids_mostrados = config.get("vertidas_mostradas")
    if args.mostrar_vertidas is not None:
        ids_mostrados = [int(item.strip()) for item in args.mostrar_vertidas.split(",") if item.strip()]

    met_e09 = metricas_tmax(resultado)
    puntos_e09 = seleccionar_puntos_representativos(grilla, resultado)
    ijk_pico = next(p.ijk for p in puntos_e09 if p.categoria == "pico_global")
    print("\nE09 · puntos representativos:")
    for punto in puntos_e09:
        print(f"  {punto.etiqueta} {punto.categoria}: ijk={punto.ijk}, "
              f"xyz={tuple(round(v, 3) for v in punto.xyz_m)} m, "
              f"vertida={punto.vertida_id}")
    print(f"E09 · Tmax global={met_e09.tmax_pico_c:.6f} °C a "
          f"t={met_e09.t_pico_h:.3f} h; Tmax final={met_e09.tmax_final_c:.6f} °C")

    vertidas_vistas = [
        v for v in vertidas
        if ids_mostrados is None or v.vertida_id in ids_mostrados
    ]
    figuras_e09 = [
        dibujar_geometria_ifc(
            vertidas, evidencia / "geometria_bim_ifc_3d.png",
            vertidas_mostradas=ids_mostrados,
        ),
        dibujar_geometria_ifc_explotada(
            vertidas_vistas,
            evidencia / "geometria_bim_ifc_3d_explotada.png",
        ),
        dibujar_cortes_termicos(
            grilla, resultado, ijk_pico,
            [met_e09.t_pico_s, resultado.tiempos_s[-1]],
            evidencia / "cortes_termicos_3d.png",
        ),
        dibujar_puntos_representativos(
            grilla, resultado, puntos_e09,
            evidencia / "temperaturas_puntos_representativos.png",
        ),
        dibujar_historial(resultado, evidencia / "tmax_historial.png"),
    ]
    perfil_vertical, columna_perfil = calcular_perfil_vertical_nucleo_superficie(
        grilla, resultado)
    guardar_perfil_vertical_e09(
        evidencia / "perfil_vertical_e09.csv", perfil_vertical)
    delta_max = max(perfil_vertical, key=lambda fila: fila[7])
    print(f"E09 · columna vertical ij={columna_perfil}; "
          f"DeltaT nucleo-superficie maximo={delta_max[7]:.6f} °C "
          f"a t={delta_max[0]:.3f} h")
    guardar_metricas_e09(
        evidencia / "metricas_e09.csv",
        grilla=grilla, resultado=resultado, puntos=puntos_e09,
        hidratacion=hidratacion, config=config, dx=dx,
        area_exterior_m2=area_exterior_m2, interfaces=interfaces,
        perfil_vertical=(perfil_vertical, columna_perfil),
    )
    guardar_puntos_e09(
        evidencia / "puntos_representativos.csv", resultado, puntos_e09,
        ruta_series=evidencia / "puntos_representativos_series.csv",
    )
    print(f"Evidencia E09 guardada en: {evidencia}")

    if modo_e09:
        print("Modo E09: solo se generaron las cinco figuras seleccionadas "
              "y sus CSV; GIF y figuras generales omitidos.")
        print(f"\nSalidas E09 en: {salida}")
        if args.mostrar:
            plt.show()
        else:
            for figura in figuras_e09:
                plt.close(figura)
        return 0

    if not args.sin_gif:
        cantidad_cuadros = crear_gif_evolucion_3d(
            grilla, resultado, tiempos_gif_s, salida / "evolucion_temperatura_3d.gif",
        )
        print(f"GIF guardado en: {salida / 'evolucion_temperatura_3d.gif'} "
              f"({cantidad_cuadros} cuadros, cada {intervalo_gif_h:g} h)")
        cantidad_cuadros_control = crear_gif_control_termico_3d(
            grilla, resultado, control,
            limite_temperatura_c, limite_diferencia_c,
            salida / "control_termico_evolucion_3d.gif",
        )
        print(f"GIF de control térmico guardado en: "
              f"{salida / 'control_termico_evolucion_3d.gif'} "
              f"({cantidad_cuadros_control} cuadros, cada {intervalo_control_h:g} h)")
    else:
        print("GIF omitidos por --sin-gif")

    figs = [
        dibujar_geometria(grilla, ids_mostrados, salida / "geometria_voxelizada_3d.png"),
        dibujar_temperatura(grilla, resultado, salida / "temperatura_final_3d.png"),
        dibujar_control_termico(
            control, salida / "control_termico.png",
            limite_temperatura_c, limite_diferencia_c,
        ),
        dibujar_evolucion_temperatura(
            grilla, resultado, tiempos_evolucion_s,
            salida / "evolucion_temperatura_3d.png",
        ),
        dibujar_historial(resultado, salida / "historial_termico.png"),
    ]
    figs.extend(figuras_e09)
    print(f"\nSalidas guardadas en: {salida}")
    if args.mostrar or config.get("mostrar_graficos", False):
        plt.show()
    else:
        for figura in figs:
            plt.close(figura)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
