"""Demo corta E11: clima desde un snapshot, simulación y exportación 3D (VTK).

Ejecutar con ▶️ (o ``python scripts/demo_e11_clima_real.py``). Tarda unos segundos.

Pasos que hace (y archivos que deja en ``results/e11_demo/``):
  1. Muestra los metadatos del snapshot (``clima.snapshot.carpeta`` de configs/e11_clima_real.jsonc).
  2. Grafica la temperatura del aire y el intervalo cubierto -> ``aire_snapshot.png``.
  3. Corre la simulación corta dentro del intervalo cubierto, con h constante (config E10).
  4. Muestra en qué función se lee la serie y dónde se evalúa T_aire(t) para el flujo Robin.
  5. Exporta estados 3D -> ``vtk/estado.pvd`` (+ ``vtk/marcas.pvd``) y deja el guion del lab.

Mientras las coordenadas oficiales (sitios.yaml) sigan pendientes el snapshot es SINTÉTICO y la
demo lo dice de forma explícita: la corrida con clima real sigue pendiente.
"""
from __future__ import annotations

import inspect
import json
import math
from pathlib import Path
import sys

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import demo_pipeline_ifc_3d as pipe
from demo_e10_entorno import correr
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal import diffusion_nd
from mcplan.thermal.seguimiento import SeguimientoPuntos, proponer_puntos, resolver_puntos
from mcplan.visualization.vtk_export import (
    CODIGO_CARA_RECONSTRUIDA, CODIGO_NUCLEO, CODIGO_VOXEL_PROXIMO, Marca, exportar_estados)
from mcplan.weather import serie as modulo_serie

CONFIGS = [RAIZ / "configs" / "semana3_fundacion.jsonc",
           RAIZ / "configs" / "e10_demo.jsonc",
           RAIZ / "configs" / "e11_clima_real.jsonc"]
SALIDA = RAIZ / "results" / "e11_demo"


def cargar_config() -> dict:
    config = pipe.cargar_configuracion(CONFIGS[0])
    for ruta in CONFIGS[1:]:
        config.update(json.loads(pipe.quitar_comentarios_jsonc(ruta.read_text(encoding="utf-8"))))
    return config


def ubicacion(objeto, texto: str | None = None) -> str:
    """``archivo:linea`` de una función o de la primera línea que contiene ``texto``."""
    lineas, inicio = inspect.getsourcelines(objeto)
    archivo = Path(inspect.getsourcefile(objeto)).relative_to(RAIZ)
    if texto:
        for n, linea in enumerate(lineas):
            if texto in linea:
                return f"{archivo}:{inicio + n}"
    return f"{archivo}:{inicio}"


def metadatos_snapshot(carpeta: Path) -> dict:
    return json.loads((carpeta / "metadata.json").read_text(encoding="utf-8"))


def graficar_aire(aire, config, final_s, ruta):
    serie = aire.serie
    inicio = aire.cobertura.solicitado[0]
    fin_sim = aire.cobertura.solicitado[1]
    fig, eje = plt.subplots(figsize=(11, 4.6))
    eje.plot(serie.tiempos_utc, serie.temperatura_c, color="#1f77b4", linewidth=1.4,
             label="Temperatura del aire del snapshot")
    eje.axvspan(inicio, fin_sim, color="#ffdd99", alpha=0.45,
                label=f"Ventana simulada ({final_s / 3600:g} h desde el inicio de obra)")
    eje.axvline(serie.intervalo[0], color="#555555", linestyle=":")
    eje.axvline(serie.intervalo[1], color="#555555", linestyle=":", label="Intervalo cubierto por el snapshot")
    tipo = aire.descripcion.get("tipo", "sintetica" if aire.es_sintetica else "pronostico")
    aviso = "  —  SERIE SINTÉTICA, NO es un dato real" if aire.es_sintetica else ""
    eje.set_title(f"Aire (snapshot {aire.descripcion.get('snapshot', '?')}, tipo: {tipo}){aviso}")
    eje.set_ylabel("T aire [°C]")
    eje.set_xlabel("Fecha y hora (UTC)")
    eje.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m %H:%M"))
    eje.grid(True, alpha=0.3)
    eje.legend(fontsize=8, loc="lower right")
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(ruta, dpi=140)
    plt.close(fig)


def main() -> int:
    SALIDA.mkdir(parents=True, exist_ok=True)
    config = cargar_config()
    cfg_snapshot = config["clima"]["snapshot"]
    carpeta = RAIZ / cfg_snapshot["carpeta"]
    print("== 1. Archivo y metadatos del snapshot ==")
    print(f"Carpeta: {carpeta}")
    print("Archivos:", ", ".join(sorted(p.name for p in carpeta.iterdir())))
    meta = metadatos_snapshot(carpeta)
    for clave in ("series_kind", "source", "model", "site", "latitude", "longitude", "coordinates_source",
                  "units", "issued_at", "retrieved_at", "valid_from", "valid_to", "original_sha256"):
        print(f"  {clave}: {meta.get(clave)}")
    original = carpeta / meta.get("original_file", "original.json")
    try:
        crudo = json.loads(original.read_text(encoding="utf-8"))
        if "elevation" in crudo:
            print(f"  punto de grilla devuelto por el proveedor: lat {crudo.get('latitude')}, lon {crudo.get('longitude')}, "
                  f"elevación {crudo.get('elevation')} m (el pedido fue lat {meta.get('latitude')}, lon {meta.get('longitude')})")
    except (OSError, ValueError):
        pass
    sintetico = bool(meta.get("synthetic"))
    if sintetico:
        print("  >>> SNAPSHOT SINTÉTICO: la corrida con clima real sigue PENDIENTE (faltan las coordenadas "
              "oficiales de sitios.yaml).")

    vertidas = pipe.aplicar_overrides(extraer_vertidas((RAIZ / config["ifc"]).resolve()), config)
    grilla = voxelizar_vertidas(vertidas, float(config["paso_voxel_m"]))
    puntos = resolver_puntos(grilla, proponer_puntos(grilla, "x+"))
    seguimiento = SeguimientoPuntos(grilla, puntos)
    final_s = float(config["tiempo_final_h"]) * 3600.0

    print("\n== 2. Corrida corta con la serie de aire del snapshot ==")
    resultado, bordes, aire = correr(config, grilla, vertidas, observadores=(seguimiento,))
    print(f"Cobertura: {aire.cobertura.resumen()}")
    graficar_aire(aire, config, final_s, SALIDA / "aire_snapshot.png")
    print(f"Curva de aire: {SALIDA / 'aire_snapshot.png'}")

    h = bordes.condicion_aire.coeficiente_h_w_m2_k
    print(f"\n== 3. Dónde se lee la serie y dónde se evalúa T_aire(t) (h constante = {h:g} W/m²K) ==")
    print(f"  carga/valida el snapshot : {ubicacion(pipe.construir_aire)} (loader.construir_aire)")
    print(f"  crea T_aire(t) [s desde el inicio de obra] : {ubicacion(modulo_serie.SerieTemperatura.a_funcion, 'def temperatura')}"
          " (interpolación lineal; fuera de la serie lanza error)")
    print(f"  se enchufa al Robin      : {ubicacion(pipe.crear_bordes_entorno, 'robin(h, aire.funcion)')} "
          "(robin(h, aire.funcion))")
    print(f"  se evalúa en cada paso   : {ubicacion(diffusion_nd.DifusionCartesiana._flujo_borde_w, 'cond.temperatura(tiempo_s)')} "
          "(cond.temperatura(tiempo_s) -> T_aire(t))")
    print("  Ecuación (Robin por cara expuesta):  Phi = G · (T_aire(t) − T_celda),   "
          "G = A / (d/(2k) + R''_protección + 1/h)   [W]")
    t_ejemplo = 36.0 * 3600.0
    punto = next(p for p in puntos if p.vertida_id == 2)
    serie_p = seguimiento.serie(2)
    n = int(np.argmin(np.abs(serie_p["tiempos_s"] - t_ejemplo)))
    t_aire = aire.funcion(float(serie_p["tiempos_s"][n]))
    t_celda = float(serie_p["celda_proxima_c"][n])
    k = float(config["materiales"]["M1"]["conductividad_w_m_k"])
    d = grilla.paso_m["xyz".index(punto.cara[0])]
    r_extra = float(bordes.resistencias_en(float(serie_p["tiempos_s"][n])).get(punto.cara, np.zeros(grilla.forma))
                    [punto.indice_voxel_proximo])
    q = (t_aire - t_celda) / (d / (2 * k) + r_extra + 1.0 / h)
    print(f"  Ejemplo V2, cara {punto.cara}, t = {serie_p['tiempos_s'][n] / 3600:.1f} h: T_aire = {t_aire:.2f} °C, "
          f"T_celda = {t_celda:.2f} °C, R'' = {r_extra:.3f} m²K/W -> q = {q:.1f} W/m² (+ entra)")

    print("\n== 4. Exportación 3D (VTK) ==")

    def marcas(n_estado, t_s):
        m = int(np.argmin(np.abs(np.asarray(seguimiento.tiempos_s) - t_s)))
        lista = []
        for p in puntos:
            s = seguimiento.serie(p.vertida_id)
            if not np.isfinite(s["interior_c"][m]):
                continue                                         # vertida aún no colocada: sin marca
            lista.append(Marca(p.vertida_id, CODIGO_NUCLEO, p.centro_interior_m, float(s["interior_c"][m])))
            lista.append(Marca(p.vertida_id, CODIGO_VOXEL_PROXIMO, p.centro_voxel_proximo_m,
                               float(s["celda_proxima_c"][m])))
            if np.isfinite(s["cara_estimada_c"][m]):
                lista.append(Marca(p.vertida_id, CODIGO_CARA_RECONSTRUIDA, p.centro_cara_m,
                                   float(s["cara_estimada_c"][m])))
        return lista

    cfg_export = config["exportacion_3d"]
    resumen = exportar_estados(RAIZ / cfg_export["carpeta"], grilla, resultado, cfg_export["estados_h"],
                               marcas_por_estado=marcas, prefijo="estado", solo_activas=True)
    for estado in resumen["estados"]:
        print(f"  t = {estado['tiempo_h']:6.1f} h: {estado['archivo']} ({estado['celdas']} celdas activas)")
    print(f"  Colección temporal: {resumen['pvd']}\n  Marcas             : {resumen['pvd_marcas']}")
    print(f"  Datos de celda: temperatura_C [°C], vertida_id, mascara_activa, n_caras_expuestas; "
          f"coordenadas en m; tiempo en FieldData (TimeValue, tiempo_h [h])")
    print("  Marcas: codigo_marca 1 = núcleo, 2 = voxel próximo a la cara (T del voxel), "
          "3 = cara exterior con T RECONSTRUIDA (distinta de la del voxel)")

    estado_aire = "SINTÉTICA (corrida real pendiente)" if sintetico else f"REAL ({meta.get('series_kind')}, {meta.get('source')})"
    resumen_txt = (f"Demo E11 · serie de aire {estado_aire}\nSnapshot: {carpeta.name} · tipo {meta.get('series_kind')}\n"
                   f"{aire.cobertura.resumen()}\nTmax pico: {resultado.temperaturas_maximas_c.max():.3f} °C\n"
                   f"h = {h:g} W/m²K; retiros: {[(round(t / 3600, 1), n) for t, n, *_ in resultado.retiros]}\n")
    (SALIDA / "resumen_e11.txt").write_text(resumen_txt, encoding="utf-8")

    print("\n== 5. Guion del lab (abrir en este orden) ==")
    print(f"  1. Archivo y metadatos : {carpeta / 'metadata.json'}  (y original.json, forecast.json)")
    print(f"  2. Curva de aire       : {SALIDA / 'aire_snapshot.png'}")
    print(f"  3. Ejecutar            : scripts/demo_e11_clima_real.py con ▶️ (esta corrida)")
    print(f"  4. Resultado 3D        : python scripts/ver_resultado_3d.py   (PyVista)   "
          f"o ParaView: abrir {SALIDA / 'vtk' / 'estado.pvd'} y {SALIDA / 'vtk' / 'marcas.pvd'}")
    print("  5. Corte               : en el visor, mover el plano de corte (slice) o activar --clip; "
          "señalar el núcleo (esfera roja) y la cara exterior reconstruida (esfera azul)")
    if sintetico:
        print("\nRECORDATORIO: esta corrida usa una serie SINTÉTICA. Con las coordenadas oficiales: "
              "scripts/descargar_snapshot_e11.py y actualizar clima.snapshot.carpeta.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
