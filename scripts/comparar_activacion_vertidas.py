"""Piloto breve: activación instantánea vs colocación progresiva por vertida.

Ejecutar desde la raíz: ``python scripts/comparar_activacion_vertidas.py``.
Usa el plan propuesto de 6 h, IFC, M1_SS, clima GEFS de referencia y el mismo
modelo de bordes E10. La progresión reparte uniformemente el volumen de cada
vertida en grupos horarios; no representa un frente real de descarga.
"""
from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import demo_e10_entorno as demo
import demo_pipeline_ifc_3d as pipe
import evaluar_dossier_final as dossier_tools
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision
from mcplan.thermal.hidratacion import ModeloFuenteTabulada


def tiempos_progresivos(grilla, vertidas, plan: PlanConstruccion) -> np.ndarray:
    """Activa grupos horarios espacialmente distribuidos dentro de cada vertida."""
    tiempos = np.full(grilla.forma, np.inf, dtype=float)
    for vertida in vertidas:
        celdas = np.flatnonzero(grilla.vertida_id == vertida.vertida_id)
        inicio_h = plan.activaciones_h[vertida.vertida_id]
        duracion_h = plan.duraciones_h[vertida.vertida_id]
        if duracion_h <= 0.0:
            tiempos.reshape(-1)[celdas] = inicio_h * 3600.0
            continue

        intervalos_h = min(6, max(1, int(np.ceil(duracion_h))))
        momentos = np.linspace(inicio_h, inicio_h + duracion_h, intervalos_h + 1)
        # Semilla fija para que la distribución espacial sea reproducible y no
        # concentre una capa o esquina completa en un solo grupo de activación.
        rng = np.random.default_rng(20261116 + int(vertida.vertida_id))
        celdas_ordenadas = rng.permutation(celdas)
        grupos = np.array_split(celdas_ordenadas, len(momentos))
        for hora, grupo in zip(momentos, grupos):
            tiempos.reshape(-1)[grupo] = float(hora) * 3600.0
    return tiempos


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ruta_cfg = RAIZ / "configs/evaluacion_final.jsonc"
    config, dossier = dossier_tools.cargar_configuracion_dossier(ruta_cfg)
    config["tiempo_final_h"] = 336.0
    config["paso_voxel_m"] = 0.5  # Malla gruesa para una comparación piloto rápida.
    config["paso_tiempo_max_s"] = 600.0

    plan_cfg = next(p for p in dossier["planes"] if p["id"] == "propuesto_6h")
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan = dossier_tools.cargar_plan(plan_cfg, originales)
    vertidas = plan.aplicar(originales)
    grilla = voxelizar_vertidas(vertidas, config["paso_voxel_m"])

    curva_cfg = next(c for c in dossier["curvas"] if c["id"] == "M1_SS")
    hidratacion = ModeloFuenteTabulada.desde_source_terms(
        RAIZ / curva_cfg["archivo_source_terms"], curva_cfg["columna"],
        cola_cero=bool(curva_cfg.get("cola_cero", True)))
    clima_cfg = next(c for c in dossier["climas"]
                     if c["id"] == "GEFS_20261116_representativo")
    clima_overlay = dossier_tools.leer_jsonc(RAIZ / clima_cfg["config_jsonc"])
    config["clima"] = dict(clima_overlay["clima"])
    config["clima"]["inicio_obra_utc"] = dossier["inicio_obra_utc"]

    escenarios = [
        ("Instantánea", None),
        ("Progresiva (grupos horarios)", tiempos_progresivos(grilla, vertidas, plan)),
    ]
    filas = []
    for etiqueta, activaciones in escenarios:
        observador = SeguimientoMetricasDecision(
            grilla, float(dossier.get("profundidad_s25_m", 0.025)),
            tuple(dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"])))
        resultado, bordes, aire = demo.correr(
            config, grilla, vertidas, observadores=(observador,),
            modelo_hidratacion=hidratacion,
            tiempos_activacion_s_por_voxel=activaciones)
        metricas = observador.resumen()
        filas.append({
            "escenario": etiqueta,
            "dx_m": grilla.paso_m[0],
            "temperatura_maxima_c": metricas["M_T_c"],
            "hora_temperatura_maxima": metricas["t_M_T_h"],
            "diferencia_maxima_interior_superficie_c": metricas["M_delta_c"],
            "hora_diferencia_maxima": metricas["t_M_delta_h"],
            "energia_balance_residuo_j": resultado.error_balance_j,
            "clima": aire.descripcion.get("snapshot") if aire else "constante",
        })
        print(f"{etiqueta}: Tmax={metricas['M_T_c']:.3f} °C "
              f"(t={metricas['t_M_T_h']:.1f} h); "
              f"Diferencia maxima interior-superficie={metricas['M_delta_c']:.3f} °C "
              f"(t={metricas['t_M_delta_h']:.1f} h); "
              f"balance={resultado.error_balance_j:.3e} J")

    base, progresiva = filas
    print("\nCambio progresiva − instantánea: "
          f"Tmax={progresiva['temperatura_maxima_c'] - base['temperatura_maxima_c']:+.3f} °C; "
          f"diferencia interior-superficie={progresiva['diferencia_maxima_interior_superficie_c'] - base['diferencia_maxima_interior_superficie_c']:+.3f} °C")
    destino = RAIZ / "results/dossier_final/comparacion_activacion_vertidas.csv"
    destino.parent.mkdir(parents=True, exist_ok=True)
    with destino.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=list(filas[0]))
        escritor.writeheader()
        escritor.writerows(filas)
    print(f"Malla: {grilla.forma}, Δx={grilla.paso_m[0]:g} m; horizonte: 336 h; "
          f"curva: {curva_cfg['id']}; clima: {clima_cfg['id']}.")
    print(f"Comparación guardada en: {destino}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
