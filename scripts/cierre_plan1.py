"""Continuación reproducible del caso 4, con un resultado atómico por escenario."""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import multiprocessing
import platform
from pathlib import Path
import resource
import subprocess
import time

import bateria7_casos as bateria
import demo_e10_entorno as demo
import numpy as np

from defensa_comprobaciones import tiempos_frente
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import SeguimientoMetricasDecision

ROOT = bateria.RAIZ
CONFIG = ROOT / "configs/defensa/cierre_plan1.json"
SPEC = json.loads(CONFIG.read_text())
SALIDA = ROOT / "results/cierre_plan1"
TABLAS = ROOT / "reports/Defensa/cierre_plan1"


class Seguimiento(SeguimientoMetricasDecision):
    """Conserva la razón de estabilidad, además de los extremos S25."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.razon_max = 0.0
        self.paso_max_s = 0.0
        self.n_pasos = 0

    def registrar_paso(self, tiempo_s, dt_s, solver, temperatura_c, fuente_w_m3):
        super().registrar_paso(tiempo_s, dt_s, solver, temperatura_c, fuente_w_m3)
        self.razon_max = max(self.razon_max, dt_s / solver.paso_maximo_estable_s())
        self.paso_max_s = max(self.paso_max_s, dt_s)
        self.n_pasos += 1


def correr(tarea):
    variante, escenario = tarea
    g = bateria._contexto()
    ev, dossier = g["ev"], g["dossier"]
    parametros = SPEC["variantes"][variante]
    _, (curva, fuente), clima, peso = g["matriz"][escenario]
    config = copy.deepcopy(g["config"])
    config.update(copy.deepcopy(g["overlays"][clima["id"]]))
    etiqueta = Path(clima["config_jsonc"]).name
    ruta_clima = ROOT / parametros.get("climas_extendidos", SPEC["climas_extendidos"]) / etiqueta
    config["clima"] = ev.leer_jsonc(ruta_clima)["clima"]
    horizonte_h = parametros.get("horizonte_h", SPEC["horizonte_h"])
    if "extension_clima" in parametros:
        config["clima"]["snapshot"]["extension"] = parametros["extension_clima"]
    config.update(tiempo_final_h=horizonte_h, paso_voxel_m=parametros["dx_m"],
                  paso_tiempo_max_s=parametros["dt_s"])
    pcfg = copy.deepcopy(g["planes_cfg"][SPEC["plan_id"]])
    if parametros.get("mantener_pur"):
        for protecciones in pcfg["entorno"]["protecciones_por_vertida"].values():
            manta = next(p for p in protecciones if p["nombre"] == "manta_PUR")
            manta.pop("retiro_h")
            manta.pop("retiro_desde")
            if "espesor_pur_mm" in parametros:
                manta["espesor_m"] = parametros["espesor_pur_mm"] / 1000
    if "retiro_capas_h" in parametros:
        retiros = parametros["retiro_capas_h"]
        for protecciones in pcfg["entorno"]["protecciones_por_vertida"].values():
            manta = next(p for p in protecciones if p["nombre"] == "manta_PUR")
            protecciones.remove(manta)
            for capa, retiro in enumerate(retiros, 1):
                protecciones.append(dict(manta, nombre=f"manta_PUR_{capa}",
                                         espesor_m=manta["espesor_m"] / len(retiros),
                                         retiro_h=retiro, retiro_desde="inicio_obra"))
    config["entorno"] = pcfg["entorno"]
    if "caudal_m3_h" in parametros:
        inicio = 0.0
        for v in sorted(g["originales"], key=lambda v: v.vertida_id):
            i = str(v.vertida_id)
            duracion = v.volumen_m3 / parametros["caudal_m3_h"]
            pcfg["activaciones_h"][i], pcfg["duraciones_h"][i] = inicio, duracion
            inicio += duracion
    plan = ev.cargar_plan(pcfg, g["originales"])
    entradas = {"config": config, "plan": pcfg, "escenario": escenario, "peso": peso,
                "curva": curva, "variante": variante}
    digest = hashlib.sha256(json.dumps(entradas, sort_keys=True).encode())
    archivos = [Path(__file__), ROOT / SPEC["bateria"], ROOT / config["ifc"],
                ROOT / curva["archivo_source_terms"], ruta_clima]
    carpeta_clima = ROOT / config["clima"]["snapshot"]["carpeta"]
    archivos.extend(carpeta_clima / nombre for nombre in ("metadata.json", "forecast.json"))
    archivos.extend(sorted((ROOT / "src/mcplan").rglob("*.py")))
    archivos.extend(ROOT / "scripts" / nombre for nombre in (
        "demo_e10_entorno.py", "demo_pipeline_ifc_3d.py", "defensa_comprobaciones.py",
        "evaluar_dossier_final.py", "bateria7_casos.py"))
    for archivo in archivos:
        digest.update(hashlib.sha256(archivo.read_bytes()).digest())
    firma = digest.hexdigest()
    destino = SALIDA / variante / f"{escenario}.json.gz"
    if destino.exists():
        with gzip.open(destino, "rt") as f:
            previo = json.load(f)
        if previo["firma"] == firma:
            return previo["fila"]

    inicio = time.perf_counter()
    vertidas = plan.aplicar(g["originales"])
    clave = (variante, parametros["dx_m"])
    if clave not in g["grillas"]:
        g["grillas"][clave] = voxelizar_vertidas(vertidas, parametros["dx_m"])
    grilla = g["grillas"][clave]
    tiempos = None
    if "caudal_m3_h" in parametros:
        tiempos = tiempos_frente(grilla, vertidas, {
            "activaciones_h": plan.activaciones_h, "duraciones_h": plan.duraciones_h})
    obs = Seguimiento(grilla, dossier["profundidad_s25_m"], tuple(dossier["caras_s25"]),
                       guardar_historial=True)
    res, _, aire = demo.correr(config, grilla, vertidas, observadores=(obs,),
                              modelo_hidratacion=fuente, guardar_campos=False,
                              tiempos_activacion_s_por_voxel=tiempos)
    if aire is None or aire.es_sintetica or obs.razon_max > 1 or not obs.n_pasos:
        raise ValueError("Corrida sin trayectoria GEFS o con verificación de estabilidad fallida")
    metricas = obs.resumen()
    original = [x for x in obs.historial if x["tiempo_h"] <= 336]
    falla_t = metricas["M_T_c"] > 70
    falla_delta = metricas["M_delta_c"] > 20
    fila = dict(variante=variante, escenario=escenario, peso=peso,
                D_h=plan.duracion_obra_h(0), dx_m=parametros["dx_m"],
                dt_limite_s=parametros["dt_s"], horizonte_h=horizonte_h,
                MT_336_c=max(x["temperatura_maxima_c"] for x in original),
                Mdelta_336_c=max(x["diferencia_maxima_S25_c"] for x in original),
                **metricas, falla_T=falla_t, falla_delta=falla_delta,
                falla_conjunta=falla_t or falla_delta, pasos=obs.n_pasos,
                dt_aplicado_max_s=obs.paso_max_s, dt_sobre_estable_max=obs.razon_max,
                energia_fuente_J=res.energia_fuente_j, energia_activacion_J=res.energia_activacion_j,
                energia_frontera_J=res.energia_frontera_j, energia_final_J=res.energia_final_j,
                balance_residual_J=res.error_balance_j,
                balance_relativo=abs(res.error_balance_j) / max(
                    1, abs(res.energia_fuente_j) + abs(res.energia_activacion_j)
                    + abs(res.energia_frontera_j)),
                segundos=time.perf_counter()-inicio)
    destino.parent.mkdir(parents=True, exist_ok=True)
    temporal = destino.with_suffix(".tmp")
    datos = dict(firma=firma, fila=fila, entradas=entradas, retiros=res.retiros,
                 serie=obs.historial, rss_max_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 python=platform.python_version(), numpy=np.__version__, plataforma=platform.platform())
    with gzip.open(temporal, "wt") as f:
        json.dump(datos, f, allow_nan=False)
    temporal.replace(destino)
    return fila


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variante", choices=SPEC["variantes"])
    parser.add_argument("--escenario", action="append")
    parser.add_argument("--procs", type=int, default=1)
    args = parser.parse_args()
    matriz = bateria._contexto()["matriz"]
    escenarios = args.escenario or list(matriz)
    tareas = [(args.variante, e) for e in escenarios]
    if args.procs == 1:
        for tarea in tareas:
            print(json.dumps(correr(tarea)), flush=True)
    else:
        with multiprocessing.get_context("fork").Pool(args.procs) as pool:
            for fila in pool.imap_unordered(correr, tareas):
                print(json.dumps(fila), flush=True)


if __name__ == "__main__":
    main()
