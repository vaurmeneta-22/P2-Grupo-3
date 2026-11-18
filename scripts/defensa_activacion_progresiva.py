"""Sensibilidad a la hipótesis de activación instantánea de cada vertida.

Compara, para un plan, tres tratamientos de la colocación física (D no cambia):
  inicio      : toda la vertida se activa al INICIO de su colocación (hipótesis de la búsqueda),
  progresiva  : los voxeles se reparten en grupos horarios durante la colocación,
  fin         : toda la vertida se activa al FIN de su colocación.

    python scripts/defensa_activacion_progresiva.py f_P2_L2m_U63 M4_SS fria
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "scripts"))

import defensa_busqueda as db  # noqa: E402
import demo_e10_entorno as demo  # noqa: E402
from comparar_activacion_vertidas import tiempos_progresivos  # noqa: E402
from mcplan.geometry.ifc import extraer_vertidas  # noqa: E402
from mcplan.geometry.voxelizacion import voxelizar_vertidas  # noqa: E402
from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision  # noqa: E402
from mcplan.thermal.hidratacion import ModeloFuenteTabulada  # noqa: E402


def correr(plan_id, curva, clima, modo, horizonte_h=1344.0, archivo_planes="configs/defensa/planes_final.json"):
    spec = next(p for p in db._leer_jsonc(RAIZ / archivo_planes)["planes"] if p["id"] == plan_id)
    hor = db.construir_plan(spec)
    config = {}
    for capa in ("configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"):
        config.update(db._leer_jsonc(RAIZ / capa))
    config.pop("sobrescribir_vertidas", None)
    config.update({"tiempo_final_h": horizonte_h, "paso_voxel_m": 0.25, "paso_tiempo_max_s": 600.0,
                   "entorno": db.construir_entorno(spec, hor)})
    config["clima"] = dict(db._leer_jsonc(RAIZ / db.CLIMAS[clima][0])["clima"])
    config["clima"]["inicio_obra_utc"] = "2026-11-19T13:00:00Z"
    if horizonte_h > 672.0:
        config["clima"]["snapshot"] = dict(config["clima"]["snapshot"], extension="ciclo_diario")
    config["sitio"] = {"grupo": 3, "latitud": -31.1115, "longitud": -71.6147}
    act, dur = dict(hor["activaciones_h"]), dict(hor["duraciones_h"])
    if modo == "fin":
        act_efectiva = {i: act[i] + dur[i] for i in act}
    else:
        act_efectiva = act
    tcol = spec.get("temperatura_colocacion_c", 18.0)
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan_fisico = PlanConstruccion(plan_id, act, dur, {i: float(tcol) for i in range(1, 6)}, frentes_simultaneos_max=1)
    plan_activ = PlanConstruccion(plan_id, act_efectiva, {i: 0.0 for i in dur}, {i: float(tcol) for i in range(1, 6)}, frentes_simultaneos_max=1)
    vertidas = plan_activ.aplicar(originales)
    grilla = voxelizar_vertidas(vertidas, 0.25)
    # Las protecciones se construyen con las horas efectivas de activación (instalación relativa a la colocación).
    config["entorno"] = db.construir_entorno(spec, {"activaciones_h": act_efectiva, "duraciones_h": dur})
    tiempos = tiempos_progresivos(grilla, plan_fisico.aplicar(originales), plan_fisico) if modo == "progresiva" else None
    obs = SeguimientoMetricasDecision(grilla, 0.025)
    obs14 = SeguimientoMetricasDecision(grilla, 0.025, tiempo_max_s=336 * 3600.0)
    fuente = ModeloFuenteTabulada.desde_source_terms(db.SOURCE_TERMS, curva, cola_cero=True)
    res, _, _ = demo.correr(config, grilla, vertidas, observadores=(obs, obs14), modelo_hidratacion=fuente,
                            tiempos_activacion_s_por_voxel=tiempos)
    m, m14 = obs.resumen(), obs14.resumen()
    return {"plan": plan_id, "curva": curva, "clima": clima, "activacion": modo,
            "M_T_14d_c": m14["M_T_c"], "M_delta_14d_c": m14["M_delta_c"],
            "M_T_fin_c": m["M_T_c"], "M_delta_fin_c": m["M_delta_c"], "t_M_delta_fin_h": m["t_M_delta_h"],
            "balance_error_J": float(res.error_balance_j)}


def main() -> int:
    plan_id, curva, clima = sys.argv[1:4]
    modos = sys.argv[4].split(",") if len(sys.argv) > 4 else ("inicio", "progresiva", "fin")
    filas = []
    for modo in modos:
        f = correr(plan_id, curva, clima, modo)
        print(json.dumps(f), flush=True)
        filas.append(f)
    (RAIZ / "results/defensa").mkdir(parents=True, exist_ok=True)
    ruta = RAIZ / "results/defensa" / f"sens_activacion_{plan_id}_{curva}_{clima}{'' if len(modos) == 3 else '_' + '_'.join(modos)}.csv"
    with ruta.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
