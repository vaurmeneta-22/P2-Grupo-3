"""Evidencia de balance de energía de un caso completo (independiente de las métricas térmicas).

Comprueba tres cosas para el plan/curva/clima indicados:
  1. E_fuente del solver = Σ V_voxel · ∫Qdot dt de la curva (calor total esperado de la hidratación).
  2. E_final = E_activación + E_fuente + E_frontera (cierre del balance, con E_final recalculada
     aquí a partir del campo final, no tomada del solver).
  3. Signo y orden de magnitud del calor perdido al entorno.

    python scripts/defensa_balance_energia.py f_P2_L2m_U63 M4_SS fria
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
import demo_pipeline_ifc_3d as pipe  # noqa: E402
from mcplan.geometry.ifc import extraer_vertidas  # noqa: E402
from mcplan.geometry.voxelizacion import voxelizar_vertidas  # noqa: E402
from mcplan.thermal.decision import PlanConstruccion  # noqa: E402
from mcplan.thermal.hidratacion import ModeloFuenteTabulada  # noqa: E402


def main() -> int:
    plan_id, curva, clima = sys.argv[1:4]
    horizonte = float(sys.argv[4]) if len(sys.argv) > 4 else 672.0
    archivo = sys.argv[5] if len(sys.argv) > 5 else "configs/defensa/planes_final.json"
    spec = next(p for p in db._leer_jsonc(RAIZ / archivo)["planes"] if p["id"] == plan_id)
    hor = db.construir_plan(spec)
    config = {}
    for capa in ("configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"):
        config.update(db._leer_jsonc(RAIZ / capa))
    config.pop("sobrescribir_vertidas", None)
    config.update({"tiempo_final_h": horizonte, "paso_voxel_m": 0.25, "paso_tiempo_max_s": 600.0,
                   "entorno": db.construir_entorno(spec, hor)})
    config["clima"] = dict(db._leer_jsonc(RAIZ / db.CLIMAS[clima][0])["clima"])
    config["clima"]["inicio_obra_utc"] = "2026-11-19T13:00:00Z"
    config["sitio"] = {"grupo": 3, "latitud": -31.1115, "longitud": -71.6147}
    tcol = float(spec.get("temperatura_colocacion_c", 18.0))
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan = PlanConstruccion(plan_id, hor["activaciones_h"], hor["duraciones_h"], {i: tcol for i in range(1, 6)},
                            frentes_simultaneos_max=int(spec.get("frentes_max", 1)))
    vertidas = plan.aplicar(originales)
    grilla = voxelizar_vertidas(vertidas, 0.25)
    fuente = ModeloFuenteTabulada.desde_source_terms(db.SOURCE_TERMS, curva, cola_cero=True)
    res, _, _ = demo.correr(config, grilla, vertidas, modelo_hidratacion=fuente)

    k, rho, cp = pipe.campos_materiales(grilla, vertidas, config["materiales"])
    vol = float(np.prod(grilla.paso_m))
    activa = np.asarray(res.mascaras_activas[-1], dtype=bool)
    campo = np.asarray(res.temperaturas_c[-1])
    e_final_indep = float(np.sum(rho[activa] * cp[activa] * campo[activa]) * vol)
    n_vox = int((grilla.vertida_id > 0).sum())
    energia_curva = float(np.trapezoid(fuente.qdot_w_m3, fuente.tiempos_s))   # J/m3 por vertida completa
    esperado = n_vox * vol * energia_curva
    cierre = e_final_indep - (res.energia_activacion_j + res.energia_fuente_j + res.energia_frontera_j)
    out = {"plan": plan_id, "curva": curva, "clima": clima, "horizonte_h": horizonte,
           "voxeles_activos": n_vox, "volumen_voxelizado_m3": n_vox * vol,
           "E_fuente_solver_J": res.energia_fuente_j, "E_fuente_esperada_J": esperado,
           "dif_fuente_relativa": res.energia_fuente_j / esperado - 1.0,
           "E_activacion_J": res.energia_activacion_j, "E_frontera_J": res.energia_frontera_j,
           "E_final_solver_J": res.energia_final_j, "E_final_independiente_J": e_final_indep,
           "dif_E_final_relativa": e_final_indep / res.energia_final_j - 1.0,
           "residuo_cierre_J": cierre, "residuo_cierre_relativo_a_fuente": cierre / res.energia_fuente_j,
           "error_balance_solver_J": res.error_balance_j,
           "fraccion_calor_perdido_al_entorno": -res.energia_frontera_j / res.energia_fuente_j}
    print(json.dumps(out, indent=1))
    ruta = RAIZ / "results/defensa/balance_energia.csv"
    ruta.parent.mkdir(parents=True, exist_ok=True)
    nuevo = not ruta.exists()
    with ruta.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out))
        if nuevo:
            w.writeheader()
        w.writerow(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
