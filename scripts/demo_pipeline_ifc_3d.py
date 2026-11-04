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
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from mcplan.geometry.ifc import VertidaIFC, extraer_vertidas
from mcplan.geometry.voxelizacion import (
    area_interfaces_entre_vertidas,
    caras_superficie_exterior,
    voxelizar_vertidas,
)
from mcplan.thermal.diffusion_nd import aislado, dirichlet, neumann, robin
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.vertidas import ResultadoEtapas, simular_vertidas


CONFIG_PREDETERMINADA = RAIZ / "configs" / "semana3_fundacion.jsonc"
SALIDA = RAIZ / "results" / "semana3"
COLORES_VERTIDA = {1: "#2c7fb8", 2: "#41b6c4", 3: "#a1dab4", 4: "#fed976", 5: "#f03b20"}


def argumentos():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG_PREDETERMINADA,
                        help="JSONC de entradas (default: configs/semana3_fundacion.jsonc)")
    parser.add_argument("--dx", type=float, help="Sobrescribe el paso de voxel en metros")
    parser.add_argument("--horas", type=float, help="Sobrescribe el horizonte de simulación")
    parser.add_argument("--mostrar-vertidas", type=str,
                        help="IDs separados por coma para la imagen geométrica; ejemplo: 1,2")
    parser.add_argument("--mostrar", action="store_true", help="Abre las figuras al terminar")
    return parser.parse_args()


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
    horas = resultado.tiempos_s / 3600.0
    fig, eje = plt.subplots(figsize=(8.2, 4.8))
    eje.plot(horas, resultado.temperaturas_maximas_c, marker="o", label="Máxima en hormigón activo")
    eje.plot(horas, resultado.temperaturas_promedio_c, marker="s", label="Media en hormigón activo")
    eje.set_xlabel("Tiempo desde primera vertida [h]")
    eje.set_ylabel("Temperatura [°C]")
    eje.set_title("Evolución térmica por activación de vertidas")
    eje.grid(True, alpha=0.35)
    eje.legend()
    fig.tight_layout()
    fig.savefig(ruta, dpi=170)
    return fig


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
    config = cargar_configuracion(args.config)
    ruta_ifc = (RAIZ / config["ifc"]).resolve()
    vertidas = aplicar_overrides(extraer_vertidas(ruta_ifc), config)
    print(f"IFC: {ruta_ifc}\nVertidas extraídas: {len(vertidas)}")
    for v in vertidas:
        print(f"  ID {v.vertida_id} {v.nombre}: inicio={v.inicio_h:g} h, "
              f"Tcol={v.temperatura_colocacion_c:g} °C, "
              f"dim={np.round(v.dimensiones_m, 3)} m, V={v.volumen_m3:.3f} m³")

    SALIDA.mkdir(parents=True, exist_ok=True)
    resoluciones = config.get("resoluciones_error_volumen_m", [])
    if resoluciones:
        print("\nConvergencia geométrica de volumen:")
        guardar_convergencia_volumen(vertidas, resoluciones, SALIDA / "error_volumen.csv")

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
        control, SALIDA / "control_termico.csv",
        limite_temperatura_c, limite_diferencia_c,
    )

    (SALIDA / "resumen.txt").write_text(
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
    guardar_historial(resultado, SALIDA / "historial_termico.csv")
    cantidad_cuadros = crear_gif_evolucion_3d(
        grilla, resultado, tiempos_gif_s, SALIDA / "evolucion_temperatura_3d.gif",
    )
    print(f"GIF guardado en: {SALIDA / 'evolucion_temperatura_3d.gif'} "
          f"({cantidad_cuadros} cuadros, cada {intervalo_gif_h:g} h)")
    cantidad_cuadros_control = crear_gif_control_termico_3d(
        grilla, resultado, control,
        limite_temperatura_c, limite_diferencia_c,
        SALIDA / "control_termico_evolucion_3d.gif",
    )
    print(f"GIF de control térmico guardado en: "
          f"{SALIDA / 'control_termico_evolucion_3d.gif'} "
          f"({cantidad_cuadros_control} cuadros, cada {intervalo_control_h:g} h)")

    ids_mostrados = config.get("vertidas_mostradas")
    if args.mostrar_vertidas is not None:
        ids_mostrados = [int(item.strip()) for item in args.mostrar_vertidas.split(",") if item.strip()]
    figs = [
        dibujar_geometria(grilla, ids_mostrados, SALIDA / "geometria_voxelizada_3d.png"),
        dibujar_temperatura(grilla, resultado, SALIDA / "temperatura_final_3d.png"),
        dibujar_control_termico(
            control, SALIDA / "control_termico.png",
            limite_temperatura_c, limite_diferencia_c,
        ),
        dibujar_evolucion_temperatura(
            grilla, resultado, tiempos_evolucion_s,
            SALIDA / "evolucion_temperatura_3d.png",
        ),
        dibujar_historial(resultado, SALIDA / "historial_termico.png"),
    ]
    print(f"\nSalidas guardadas en: {SALIDA}")
    if args.mostrar or config.get("mostrar_graficos", False):
        plt.show()
    else:
        for figura in figs:
            plt.close(figura)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
