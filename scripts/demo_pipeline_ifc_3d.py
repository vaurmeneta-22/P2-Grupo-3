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
    tiempos_evolucion_s = [float(hora) * 3600.0
                           for hora in config.get("tiempos_evolucion_h", [])
                           if float(hora) * 3600.0 <= final_s]
    if not tiempos_evolucion_s or tiempos_evolucion_s[0] != 0.0:
        tiempos_evolucion_s.insert(0, 0.0)
    if not np.isclose(tiempos_evolucion_s[-1], final_s):
        tiempos_evolucion_s.append(final_s)
    intervalo_gif_h = float(config.get("intervalo_gif_h", 4.0))
    if intervalo_gif_h <= 0:
        raise ValueError("intervalo_gif_h debe ser positivo")
    intervalo_gif_s = intervalo_gif_h * 3600.0
    tiempos_gif_s = np.arange(0.0, final_s + 1e-9, intervalo_gif_s).tolist()
    if not tiempos_gif_s or not np.isclose(tiempos_gif_s[-1], final_s):
        tiempos_gif_s.append(final_s)
    tiempos_guardado_s = sorted(set(tiempos_evolucion_s + tiempos_gif_s))
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

    ids_mostrados = config.get("vertidas_mostradas")
    if args.mostrar_vertidas is not None:
        ids_mostrados = [int(item.strip()) for item in args.mostrar_vertidas.split(",") if item.strip()]
    figs = [
        dibujar_geometria(grilla, ids_mostrados, SALIDA / "geometria_voxelizada_3d.png"),
        dibujar_temperatura(grilla, resultado, SALIDA / "temperatura_final_3d.png"),
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
