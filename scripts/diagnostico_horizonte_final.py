"""Compara el cierre térmico a 14 días y documenta la hipótesis de una vertida única.

Ejecuta dos casos diagnósticos con el IFC y el clima GEFS central: el plan
propuesto y todas las partes del IFC activadas juntas en t=0. Este segundo caso
es solo una comparación térmica idealizada, no una secuencia de obra ejecutable.
"""
from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import demo_e10_entorno as demo
import demo_pipeline_ifc_3d as pipe
import evaluar_dossier_final as evaluador
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision
from mcplan.thermal.hidratacion import ModeloFuenteTabulada
from mcplan.weather.loader import construir_aire


def _cargar() -> tuple[dict, dict, dict, list]:
    config, dossier = evaluador.cargar_configuracion_dossier(
        ROOT / "configs/evaluacion_final.jsonc")
    t0 = dossier["inicio_obra_utc"]
    clima_cfg = next(c for c in dossier["climas"] if c["id"].endswith("_central"))
    overlay = evaluador.leer_jsonc(ROOT / clima_cfg["config_jsonc"])
    for campo in ("clima", "sitio", "entorno", "condicion_superficie"):
        if campo in overlay:
            config[campo] = overlay[campo]
    config["clima"]["inicio_obra_utc"] = t0
    aire = construir_aire(config["clima"], tiempo_final_s=336 * 3600, raiz_repo=ROOT)
    if aire.es_sintetica:
        raise ValueError("Se requiere el snapshot GEFS real")
    originales = extraer_vertidas((ROOT / config["ifc"]).resolve())
    return config, dossier, clima_cfg, originales


def _correr(nombre: str, config: dict, dossier: dict, originales: list,
            activaciones: dict[int, float], duraciones: dict[int, float]) -> list[dict]:
    plan = PlanConstruccion(nombre, activaciones, duraciones,
        {int(v.vertida_id): float(v.temperatura_colocacion_c) for v in originales},
        frentes_simultaneos_max=1)
    vertidas = plan.aplicar(originales)
    grilla = voxelizar_vertidas(vertidas, float(dossier["paso_voxel_m"]))
    curva = dossier["curvas"][0]
    fuente = ModeloFuenteTabulada.desde_source_terms(
        ROOT / curva["archivo_source_terms"], curva["columna"],
        cola_cero=bool(curva.get("cola_cero", True)))
    observador = SeguimientoMetricasDecision(
        grilla, float(dossier["profundidad_s25_m"]), tuple(dossier["caras_s25"]))
    resultado, _, aire = demo.correr(config, grilla, vertidas,
                                     observadores=(observador,), modelo_hidratacion=fuente)
    if aire is None or aire.es_sintetica:
        raise ValueError("La corrida no usó clima GEFS")
    metricas = observador.resumen()
    print(f"{nombre}: M_T={metricas['M_T_c']:.2f} °C; "
          f"M_delta={metricas['M_delta_c']:.2f} °C; "
          f"t(M_T)={metricas['t_M_T_h']:.1f} h")
    filas = []
    tiempos = np.asarray(resultado.tiempos_s)
    muestras = {}
    for objetivo_h in (312.0, 336.0):
        indice = int(np.argmin(np.abs(tiempos - objetivo_h * 3600.0)))
        if abs(float(tiempos[indice]) - objetivo_h * 3600.0) > 1.0:
            raise ValueError(f"No existe muestra exacta para t={objetivo_h:g} h")
        campo = np.asarray(resultado.temperaturas_c[indice])
        activa = np.asarray(resultado.mascaras_activas[indice], dtype=bool)
        muestras[objetivo_h] = {}
        for vertida_id in sorted(int(v.vertida_id) for v in originales):
            mascara = activa & (grilla.vertida_id == vertida_id)
            muestras[objetivo_h][vertida_id] = float(np.max(campo[mascara]))
    for vertida_id in sorted(int(v.vertida_id) for v in originales):
        antes, final = muestras[312.0][vertida_id], muestras[336.0][vertida_id]
        filas.append({"caso": nombre, "vertida_id": vertida_id,
            "Tmax_t312h_c": antes, "Tmax_t336h_c": final,
            "cambio_ultimas_24h_c": final - antes,
            "sigue_aumentando_mas_0_1c": final - antes > 0.1})
    return filas


def main() -> int:
    config, dossier, _, originales = _cargar()
    ids = sorted(int(v.vertida_id) for v in originales)
    propuesto = next(p for p in dossier["planes"] if p["id"] == "propuesto_6h")
    plan = evaluador.cargar_plan(propuesto, originales)
    filas = _correr("propuesto_6h", config, dossier, originales,
                    plan.activaciones_h, plan.duraciones_h)
    filas += _correr("una_vertida_idealizada_t0", config, dossier, originales,
                     {i: 0.0 for i in ids}, {i: 0.0 for i in ids})
    salida = ROOT / "results/dossier_final/diagnostico_horizonte_y_una_vertida.csv"
    with salida.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(filas[0]))
        writer.writeheader()
        writer.writerows(filas)
    print(f"\nDiagnóstico guardado: {salida}")
    print("Tmax por vertida comparado entre 312 h y 336 h (últimas 24 h).")
    print("El caso simultáneo es una idealización para comparar geometría completa, no un plan constructivo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
