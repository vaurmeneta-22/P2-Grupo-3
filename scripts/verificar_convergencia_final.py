"""Convergencia espacial/temporal del escenario de referencia del dossier.

Ejecutar desde la raíz del repo:
    python scripts/verificar_convergencia_final.py

Corre solo cuatro casos (no la matriz completa de planes × curvas) y deja el
resumen en results/dossier_final/convergencia_numerica.csv.
"""
from __future__ import annotations

import csv
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import demo_e10_entorno as demo
import evaluar_dossier_final as dossier_tools
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import SeguimientoMetricasDecision
from mcplan.weather.loader import construir_aire


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    ruta_dossier = RAIZ / "configs/evaluacion_final.jsonc"
    config, dossier = dossier_tools.cargar_configuracion_dossier(ruta_dossier)
    cfg_conv = dossier.get("convergencia_numerica", {})
    dxs = [float(v) for v in cfg_conv.get("pasos_voxel_m", [0.5, 0.25])]
    dts = [float(v) for v in cfg_conv.get("pasos_tiempo_max_s", [600.0, 300.0])]
    if len(set(dxs)) < 2 or len(set(dts)) < 2:
        raise ValueError("Se requieren al menos dos Δx y dos pasos temporales distintos")

    config["tiempo_final_h"] = float(dossier.get("tiempo_final_h", 336.0))
    config["paso_tiempo_max_s"] = max(dts)
    config["clima"] = None
    t0 = dossier.get("inicio_obra_utc")
    ruta_ifc = (RAIZ / config["ifc"]).resolve()
    originales = extraer_vertidas(ruta_ifc)

    plan_cfg = next((p for p in dossier["planes"]
                     if p["id"] == cfg_conv.get("plan_id", "propuesto_6h")), None)
    if plan_cfg is None:
        raise ValueError("No se encontró el plan configurado para la convergencia")
    plan = dossier_tools.cargar_plan(plan_cfg, originales)
    vertidas = plan.aplicar(originales)

    escenarios = dossier_tools.construir_matriz(
        dossier_tools.leer_curvas(dossier), dossier.get("climas", []), dossier)
    escenario_id_cfg = cfg_conv.get("escenario_id", "M1_SS__GEFS_20261116_representativo")
    escenario = next((item for item in escenarios if item[0] == escenario_id_cfg), None)
    if escenario is None:
        raise ValueError(f"No se encontró el escenario de referencia: {escenario_id_cfg}")
    escenario_id, (curva_cfg, hidratacion), clima_cfg, _peso = escenario

    overlay = dossier_tools.leer_jsonc((RAIZ / clima_cfg["config_jsonc"]).resolve())
    config["clima"] = dict(overlay.get("clima") or {})
    config["clima"]["inicio_obra_utc"] = t0
    aire = construir_aire(config["clima"],
                          tiempo_final_s=config["tiempo_final_h"] * 3600.0,
                          raiz_repo=RAIZ)
    if aire.es_sintetica:
        raise ValueError("La convergencia exige el snapshot climático de referencia, no una serie sintética")

    filas = []
    total = len(set(dxs)) * len(set(dts))
    indice = 0
    limite_t = float(dossier.get("limites", {}).get("temperatura_maxima_c", 70.0))
    limite_delta = float(dossier.get("limites", {}).get("diferencia_interior_superficie_c", 20.0))
    profundidad = float(dossier.get("profundidad_s25_m", 0.025))
    caras = tuple(dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"]))

    print(f"Convergencia: plan={plan.plan_id}, escenario={escenario_id}, "
          f"curva={curva_cfg['id']}, clima={clima_cfg['id']}")
    print(f"IFC: {ruta_ifc}; t0={t0}; horizonte={config['tiempo_final_h']:g} h")
    for dx in sorted(set(dxs), reverse=True):
        grilla = voxelizar_vertidas(vertidas, dx)
        for dt_max in sorted(set(dts), reverse=True):
            indice += 1
            observador = SeguimientoMetricasDecision(grilla, profundidad, caras)
            inicio = time.perf_counter()
            resultado, _bordes, aire_corrida = demo.correr(
                {**config, "paso_voxel_m": dx, "paso_tiempo_max_s": dt_max},
                grilla, vertidas, observadores=(observador,),
                modelo_hidratacion=hidratacion)
            if aire_corrida is None or aire_corrida.es_sintetica:
                raise ValueError("La corrida no usó el clima GEFS de referencia")
            met = observador.resumen()
            filas.append({
                "plan_id": plan.plan_id,
                "escenario_id": escenario_id,
                "curva_id": curva_cfg["id"],
                "clima_id": clima_cfg["id"],
                "dx_m": dx,
                "paso_tiempo_max_s": dt_max,
                "forma_grilla": "x".join(str(n) for n in grilla.forma),
                "voxeles_hormigon": int(grilla.mascara_solido.sum()),
                "M_T_c": met["M_T_c"],
                "t_M_T_h": met["t_M_T_h"],
                "M_delta_c": met["M_delta_c"],
                "t_M_delta_h": met["t_M_delta_h"],
                "clasifica_T_falla": met["M_T_c"] > limite_t,
                "clasifica_delta_falla": met["M_delta_c"] > limite_delta,
                "residuo_balance_j": resultado.error_balance_j,
                "duracion_corrida_s": time.perf_counter() - inicio,
            })
            print(f"[{indice}/{total}] dx={dx:g} m, dt_max={dt_max:g} s: "
                  f"M_T={met['M_T_c']:.3f} °C (t={met['t_M_T_h']:.1f} h), "
                  f"M_delta={met['M_delta_c']:.3f} °C (t={met['t_M_delta_h']:.1f} h), "
                  f"balance={resultado.error_balance_j:.2e} J")

    dx_fino, dt_fino = min(set(dxs)), min(set(dts))
    buscar = lambda dx, dt: next(f for f in filas if f["dx_m"] == dx and f["paso_tiempo_max_s"] == dt)
    espacial_grueso = buscar(max(set(dxs)), dt_fino)
    espacial_fino = buscar(dx_fino, dt_fino)
    temporal_grueso = buscar(dx_fino, max(set(dts)))
    temporal_fino = buscar(dx_fino, dt_fino)
    print("\nSensibilidad espacial (dt fino fijo): "
          f"ΔM_T={espacial_grueso['M_T_c'] - espacial_fino['M_T_c']:+.3f} °C; "
          f"ΔM_delta={espacial_grueso['M_delta_c'] - espacial_fino['M_delta_c']:+.3f} °C")
    print("Sensibilidad temporal (dx fino fijo): "
          f"ΔM_T={temporal_grueso['M_T_c'] - temporal_fino['M_T_c']:+.3f} °C; "
          f"ΔM_delta={temporal_grueso['M_delta_c'] - temporal_fino['M_delta_c']:+.3f} °C")
    salida = RAIZ / dossier.get("salida", "results/dossier_final") / "convergencia_numerica.csv"
    salida.parent.mkdir(parents=True, exist_ok=True)
    with salida.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=list(filas[0]))
        escritor.writeheader()
        escritor.writerows(filas)
    print(f"\nTabla guardada: {salida}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR DE CONFIGURACIÓN/DATOS: {exc}", file=sys.stderr)
        raise SystemExit(2)
