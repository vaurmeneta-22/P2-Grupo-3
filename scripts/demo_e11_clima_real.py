"""Demo corta E11: clima desde un snapshot, simulación y exportación 3D (VTK).

Ejecutar con ▶️ (o ``python scripts/demo_e11_clima_real.py``). Tarda unos segundos.

Pasos que hace (y archivos que deja en ``results/e11_demo/``):
  1. Muestra los metadatos del snapshot (``clima.snapshot.carpeta`` de configs/e11_clima_real.jsonc).
  2. Grafica la temperatura del aire y el intervalo cubierto -> ``aire_snapshot.png``.
  3. Corre la simulación corta dentro del intervalo cubierto, con h constante (config E10).
  4. Muestra en qué función se lee la serie y dónde se evalúa T_aire(t) para el flujo Robin.
  5. Figuras de evolución térmica (generación, temperaturas, flujos, diferencias) y sus CSV.
  6. Exporta estados 3D -> ``vtk/estado.pvd`` (+ ``vtk/marcas.pvd``) y deja el guion del lab.

Qué figura mirar (todas con tiempo global desde el inicio de obra, de UNA misma simulación):
  * GENERACIÓN  -> ``generacion_hidratacion.png`` (Qdot_h por vertida y potencia total) y
                   ``flujos_por_edad.png`` (mismo pulso, comparado a igual edad).
  * ACUMULACIÓN -> ``temperaturas_por_vertida.png`` y ``temperaturas_conjunto.png`` (cuánto sube la T al
                   generar más de lo que se pierde) y ``diferencias_termicas.png`` (ΔT interior-superficie).
  * PÉRDIDAS    -> ``flujos_entorno.png`` (aire protegido/expuesto, suelo, retiro de protecciones).
Opción de configuración para abrir las figuras en ventanas: ``graficos.mostrar`` (o ``--mostrar``).

Mientras las coordenadas oficiales (sitios.yaml) sigan pendientes el snapshot es SINTÉTICO y la
demo lo dice de forma explícita: la corrida con clima real sigue pendiente.
"""
from __future__ import annotations

import csv
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
from mcplan.thermal.analisis import eventos_activacion
from mcplan.thermal.extremos import SeguimientoExtremos
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.registro_flujos import RegistroEnergetico, mascara_base
from mcplan.thermal.seguimiento import SeguimientoPuntos, proponer_puntos, resolver_puntos
from mcplan.visualization import figuras_e11 as fig11
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


def graficar_aire(aire, config, final_s, ruta, *, mostrar=False):
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
    if not mostrar:
        plt.close(fig)
    return fig


def verificar_registros(reg, resultado, hidratacion, vertidas, rho, cp):
    """Comprueba que los registros coinciden con la fuente y los flujos que usó el solver."""
    inicio_s = {v.vertida_id: v.inicio_h * 3600.0 for v in vertidas}
    t = np.asarray(reg.t_s)
    max_dif = 0.0
    for i in reg.ids:
        esperado = np.where(t >= inicio_s[i], hidratacion.generacion_w_m3(t - inicio_s[i], rho, cp), np.nan)
        registrado = np.asarray(reg.qdot_w_m3[i])
        valido = np.isfinite(registrado)
        max_dif = max(max_dif, float(np.max(np.abs(registrado[valido] - esperado[valido]))) if np.any(valido) else 0.0)
    rel_fuente = abs(reg.energia_fuente_j - resultado.energia_fuente_j) / max(abs(resultado.energia_fuente_j), 1.0)
    rel_front = abs(reg.energia_frontera_j - resultado.energia_frontera_j) / max(abs(resultado.energia_frontera_j), 1.0)
    return {"dif_qdot_modelo_w_m3": max_dif, "rel_energia_fuente": rel_fuente, "rel_energia_frontera": rel_front,
            "energia_fuente_j": reg.energia_fuente_j, "energia_fuente_solver_j": resultado.energia_fuente_j,
            "energia_frontera_j": reg.energia_frontera_j, "energia_frontera_solver_j": resultado.energia_frontera_j}


def guardar_maximos(ruta, reg, seguimiento, extremos, horas_final):
    """Máximos temporales de las magnitudes de las figuras, con su hora."""
    filas = []
    for i in reg.ids:
        pico = reg.pico(reg.qdot_w_m3[i])
        if pico:
            filas.append(("generacion_hidratacion", f"Qdot_h V{i}", pico[0], pico[1], "W/m3", "máximo del pulso"))
    pico = reg.pico(reg.potencia_total_kw())
    if pico:
        filas.append(("generacion_hidratacion", "potencia total generada", pico[0], pico[1], "kW", "Σ Qdot_h · V_voxel"))
    for fila in seguimiento.resumen():
        if fila["dt_cara"]:
            r = fila["dt_cara"]
            filas.append(("diferencias_termicas", f"ΔT interior − cara reconstruida V{fila['vertida_id']}",
                          r["tiempo_h"], r["diferencia_con_signo_c"], "°C", "máx |ΔT|"))
    for fila in extremos.resumen():
        if fila["alcance"] == "global" and fila["base"] == "todos_los_pasos" and fila["valido"]:
            filas.append(("diferencias_termicas", fila["metrica"], fila["tiempo_h"], fila["abs_c"], "°C", "dominio activo"))
    for clave in reg.caras:
        pico = reg.pico(reg.cara_w_m2[clave], signo=-1.0)
        if pico:
            filas.append(("flujos_entorno", f"q V{clave[0]} {clave[1]}", pico[0], pico[1], "W/m2", "mayor salida (más negativo)"))
    if reg.celdas_suelo is not None:
        pico = reg.pico(reg.suelo_w_m2, signo=-1.0)
        if pico:
            filas.append(("flujos_entorno", "q suelo", pico[0], pico[1], "W/m2", "mayor salida (más negativo)"))
    with Path(ruta).open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["figura", "magnitud", "hora_h", "valor", "unidad", "nota"])
        for fila in filas:
            w.writerow([fila[0], fila[1], f"{fila[2]:.4f}", f"{fila[3]:.6g}", fila[4], fila[5]])
    return filas


def elegir_backend(mostrar: bool) -> None:
    """Con ``mostrar`` usa un backend gráfico (antes de crear figuras); si no, Agg (solo guarda)."""
    if not mostrar:
        return
    for backend in ("QtAgg", "TkAgg"):
        try:
            plt.switch_backend(backend)
            return
        except Exception:  # noqa: BLE001 - probar el siguiente backend disponible
            continue
    print("AVISO: no hay un backend gráfico disponible; las figuras solo se guardan.")


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
    extremos = SeguimientoExtremos(grilla)
    reg = RegistroEnergetico(grilla, caras=[(p.vertida_id, p.cara) for p in puntos],
                             celdas_suelo=mascara_base(grilla))
    final_s = float(config["tiempo_final_h"]) * 3600.0
    mostrar = bool((config.get("graficos") or {}).get("mostrar", False)) or "--mostrar" in sys.argv[1:]
    elegir_backend(mostrar)

    print("\n== 2. Corrida corta con la serie de aire del snapshot ==")
    resultado, bordes, aire = correr(config, grilla, vertidas, observadores=(seguimiento, extremos, reg))
    assert np.array_equal(reg.celdas_suelo, bordes.mascara_base), "la base del registro debe ser la del modelo de bordes"
    print(f"Cobertura: {aire.cobertura.resumen()}")
    figuras = [graficar_aire(aire, config, final_s, SALIDA / "aire_snapshot.png", mostrar=mostrar)]
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

    print("\n== 4. Evolución térmica: generación, temperaturas, flujos y diferencias ==")
    hcfg = config["hidratacion"]
    hidratacion = ModeloPulsoCosenoidal(float(hcfg["qmax_w_m3"]), float(hcfg["duracion_pulso_h"]) * 3600.0)
    material = config["materiales"][vertidas[0].tipo_hormigon]
    verif = verificar_registros(reg, resultado, hidratacion, vertidas,
                                float(material["densidad_kg_m3"]), float(material["calor_especifico_j_kg_k"]))
    print(f"  Hidratación aplicada: pulso cosenoidal Qmax = {hidratacion.qmax_w_m3:g} W/m³, duración = "
          f"{hidratacion.duracion_s / 3600:g} h por vertida (misma fuente que usa el solver).")
    print(f"  Verificación vs. solver: |Qdot registrado − modelo| máx = {verif['dif_qdot_modelo_w_m3']:.3e} W/m³; "
          f"energía generada {verif['energia_fuente_j']:.6e} J (solver {verif['energia_fuente_solver_j']:.6e} J, "
          f"dif. rel. {verif['rel_energia_fuente']:.1e}); energía por la frontera {verif['energia_frontera_j']:.6e} J "
          f"(solver {verif['energia_frontera_solver_j']:.6e} J, dif. rel. {verif['rel_energia_frontera']:.1e}).")
    assert verif["dif_qdot_modelo_w_m3"] < 1e-9 and verif["rel_energia_fuente"] < 1e-9 and verif["rel_energia_frontera"] < 1e-9
    horas_final = final_s / 3600.0
    malla_h = np.arange(0.0, horas_final + 1e-9, 0.5)
    datos = {
        "activaciones": [(h_, ids) for h_, ids in eventos_activacion(resultado) if h_ <= horas_final],
        "retiros": [(t_s / 3600.0, nombres) for t_s, nombres, *_ in resultado.retiros],
        "protecciones": list(bordes.protecciones),
        "inicio_h": {v.vertida_id: v.inicio_h for v in vertidas},
        "aire_horas": malla_h,
        "aire_c": [aire.funcion(float(t) * 3600.0) for t in malla_h],
    }
    figuras += [
        fig11.figura_hidratacion(reg, datos, horas_final, SALIDA / "generacion_hidratacion.png"),
        fig11.figura_temperaturas_por_vertida(seguimiento, datos, horas_final, SALIDA / "temperaturas_por_vertida.png"),
        fig11.figura_conjunto(extremos, datos, horas_final, SALIDA / "temperaturas_conjunto.png"),
        fig11.figura_flujos(reg, datos, horas_final, SALIDA / "flujos_entorno.png"),
        fig11.figura_por_edad(reg, datos, SALIDA / "flujos_por_edad.png"),
        fig11.figura_diferencias(seguimiento, extremos, datos, horas_final, SALIDA / "diferencias_termicas.png"),
    ]
    reg.guardar_hidratacion_csv(SALIDA / "hidratacion_serie.csv")
    reg.guardar_flujos_csv(SALIDA / "flujos_serie.csv")
    seguimiento.guardar_csv(SALIDA / "temperaturas_vertidas_serie.csv")
    seguimiento.guardar_resumen_csv(SALIDA / "temperaturas_vertidas_resumen.csv")
    extremos.guardar_evolucion_csv(SALIDA / "conjunto_extremos_serie.csv")
    with (SALIDA / "aire_serie.csv").open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(["tiempo_h", "T_aire_C"])
        for t, valor in zip(malla_h, datos["aire_c"]):
            escritor.writerow([f"{t:.3f}", f"{valor:.4f}"])
    maximos = guardar_maximos(SALIDA / "maximos_e11.csv", reg, seguimiento, extremos, horas_final)
    for figura, magnitud, hora, valor, unidad, _ in maximos:
        if figura in ("generacion_hidratacion", "diferencias_termicas") and ("total" in magnitud or "V2" in magnitud
                                                                         or "global" in magnitud or "contraste" in magnitud
                                                                         or "dT" in magnitud):
            print(f"  máx {magnitud}: {valor:.3g} {unidad} a las {hora:.1f} h")
    print("  Convenio de flujos: + entra al hormigón, − sale. W/m² = por unidad de área; W = potencia integrada.")
    print(f"  Figuras y CSV en {SALIDA} (máximos en maximos_e11.csv).")

    print("\n== 5. Exportación 3D (VTK) ==")

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

    print("\n== 6. Guion del lab (abrir en este orden) ==")
    print(f"  1. Archivo y metadatos : {carpeta / 'metadata.json'}  (y original.json, forecast.json)")
    print(f"  2. Curva de aire       : {SALIDA / 'aire_snapshot.png'}")
    print(f"  3. Ejecutar            : scripts/demo_e11_clima_real.py con ▶️ (esta corrida)")
    print(f"  4. Resultado 3D        : python scripts/ver_resultado_3d.py   (PyVista)   "
          f"o ParaView: abrir {SALIDA / 'vtk' / 'estado.pvd'} y {SALIDA / 'vtk' / 'marcas.pvd'}")
    print("  5. Corte               : en el visor, mover el plano de corte (slice) o activar --clip; "
          "señalar el núcleo (esfera roja) y la cara exterior reconstruida (esfera azul)")
    print("  Para explicar: generación -> generacion_hidratacion.png; acumulación -> temperaturas_conjunto.png y "
          "temperaturas_por_vertida.png; pérdidas -> flujos_entorno.png (y flujos_por_edad.png para comparar a igual edad).")
    if mostrar:
        plt.show()                                  # única llamada a plt.show(): abre todas las figuras
    else:
        for figura in figuras:
            plt.close(figura)
    if sintetico:
        print("\nRECORDATORIO: esta corrida usa una serie SINTÉTICA. Con las coordenadas oficiales: "
              "scripts/descargar_snapshot_e11.py y actualizar clima.snapshot.carpeta.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
