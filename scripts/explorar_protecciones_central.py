"""Compara capas de protección con M1_SS y el pronóstico GEFS central."""
from __future__ import annotations

import copy
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import demo_e10_entorno as demo
import evaluar_dossier_final as ev
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import SeguimientoMetricasDecision
from mcplan.thermal.hidratacion import ModeloFuenteTabulada
from mcplan.weather.loader import construir_aire

OUTPUT = Path.home() / "Documents/Codex/2026-10-15/h/protecciones_central_M1_SS.csv"
CARAS_MANTA = ["x-", "x+", "y-", "y+", "z+"]


def configuraciones_proteccion(entorno_base: dict, ids: list[int]) -> list[tuple[str, dict]]:
    base = copy.deepcopy(entorno_base)
    protecciones_base = base.get("protecciones_por_vertida", {})
    solo_moldaje = {
        str(i): [p for p in protecciones_base.get(str(i), [])
                 if str(p.get("nombre", "")).lower() == "moldaje"]
        for i in ids
    }

    opciones = [("actual", base)]

    superior = copy.deepcopy(base)
    superior_cfg = copy.deepcopy(protecciones_base)
    for i in ids:
        superior_cfg.setdefault(str(i), []).append({
            "nombre": "manta_superior_3cm", "caras": ["z+"],
            "espesor_m": 0.03, "conductividad_w_m_k": 0.04,
            "retiro_h": 120.0, "retiro_desde": "colocacion",
        })
    superior["protecciones_por_vertida"] = superior_cfg
    opciones.append(("manta_superior_3cm_120h", superior))

    for etiqueta, espesor, retiro_h in (
        ("manta_lateral_superior_3cm_168h", 0.03, 168.0),
        ("manta_lateral_superior_5cm_240h", 0.05, 240.0),
    ):
        entorno = copy.deepcopy(base)
        uniforme = copy.deepcopy(solo_moldaje)
        for i in ids:
            uniforme[str(i)].append({
                "nombre": "manta_curado", "caras": CARAS_MANTA,
                "espesor_m": espesor, "conductividad_w_m_k": 0.04,
                "retiro_h": retiro_h, "retiro_desde": "colocacion",
            })
        entorno["protecciones_por_vertida"] = uniforme
        opciones.append((etiqueta, entorno))

    # Sensibilidad hipotética: la cara inferior aparece como S25 mínima después
    # de proteger laterales y coronamiento. Se conserva todo lo demás y se añade
    # una capa permanente solo bajo la vertida que apoya en el suelo.
    caso_lateral_superior = opciones[-1][1]
    for espesor in (0.05, 0.10):
        entorno = copy.deepcopy(caso_lateral_superior)
        entorno["protecciones_por_vertida"]["1"].append({
            "nombre": f"aislacion_inferior_{int(espesor * 100)}cm",
            "caras": ["z-"], "espesor_m": espesor,
            "conductividad_w_m_k": 0.04,
            "retiro_h": None, "retiro_desde": "colocacion",
        })
        opciones.append((f"manta_5cm_240h_base_{int(espesor * 100)}cm_hipotetica", entorno))

    # El máximo de contraste puede desplazarse al final del horizonte, después
    # de retirar las mantas. Este caso las conserva durante los 14 días completos.
    entorno_largo = copy.deepcopy(base)
    entorno_largo["protecciones_por_vertida"] = copy.deepcopy(solo_moldaje)
    for i in ids:
        entorno_largo["protecciones_por_vertida"][str(i)].append({
            "nombre": "manta_curado_5cm_14dias", "caras": CARAS_MANTA,
            "espesor_m": 0.05, "conductividad_w_m_k": 0.04,
            "retiro_h": 360.0, "retiro_desde": "colocacion",
        })
    entorno_largo["protecciones_por_vertida"]["1"].append({
        "nombre": "aislacion_inferior_10cm", "caras": ["z-"],
        "espesor_m": 0.10, "conductividad_w_m_k": 0.04,
        "retiro_h": None, "retiro_desde": "colocacion",
    })
    opciones.append(("manta_5cm_360h_base_10cm_hipotetica", entorno_largo))

    for espesor in (0.06, 0.075):
        entorno = copy.deepcopy(base)
        uniforme = copy.deepcopy(solo_moldaje)
        for i in ids:
            uniforme[str(i)].append({
                "nombre": f"manta_curado_{espesor:g}m_14dias", "caras": CARAS_MANTA,
                "espesor_m": espesor, "conductividad_w_m_k": 0.04,
                "retiro_h": 360.0, "retiro_desde": "colocacion",
            })
        uniforme["1"].append({
            "nombre": "aislacion_inferior_10cm", "caras": ["z-"],
            "espesor_m": 0.10, "conductividad_w_m_k": 0.04,
            "retiro_h": None, "retiro_desde": "colocacion",
        })
        entorno["protecciones_por_vertida"] = uniforme
        nombre_espesor = "6cm" if espesor == 0.06 else "7_5cm"
        opciones.append((f"manta_{nombre_espesor}_360h_base_10cm_hipotetica", entorno))
    return opciones


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    config, dossier = ev.cargar_configuracion_dossier(
        ROOT / "configs/exploracion_planes_central.jsonc")
    plan_cfg = next(p for p in dossier["planes"] if p["id"] == "base_6h_18c")
    ruta_ifc = (ROOT / config["ifc"]).resolve()
    originales = extraer_vertidas(ruta_ifc)
    plan = ev.cargar_plan(plan_cfg, originales)
    vertidas = plan.aplicar(originales)
    ids = sorted(int(v.vertida_id) for v in vertidas)
    grilla = voxelizar_vertidas(vertidas, float(dossier["paso_voxel_m"]))

    climate_cfg = next(c for c in dossier["climas"] if c["id"].endswith("_central"))
    overlay = ev.leer_jsonc((ROOT / climate_cfg["config_jsonc"]).resolve())
    for campo in ("clima", "sitio", "entorno", "condicion_superficie"):
        if campo in overlay:
            config[campo] = overlay[campo]
    config["clima"]["inicio_obra_utc"] = dossier["inicio_obra_utc"]
    aire = construir_aire(config["clima"], tiempo_final_s=336 * 3600.0, raiz_repo=ROOT)
    if aire.es_sintetica:
        raise ValueError("La prueba requiere el snapshot GEFS central, no clima sintético")

    curva = dossier["curvas"][0]
    fuente = ModeloFuenteTabulada.desde_source_terms(
        ROOT / curva["archivo_source_terms"], curva["columna"],
        cola_cero=bool(curva.get("cola_cero", True)))
    opciones = configuraciones_proteccion(config.get("entorno") or {}, ids)
    if len(sys.argv) == 3 and sys.argv[1] == "--solo":
        opciones = [op for op in opciones if op[0] == sys.argv[2]]
        if not opciones:
            raise ValueError(f"Variante desconocida para --solo: {sys.argv[2]}")
    elif len(sys.argv) != 1:
        raise ValueError("Uso: explorar_protecciones_central.py [--solo ID_VARIANTE]")
    filas = []
    for indice, (etiqueta, entorno) in enumerate(opciones, start=1):
        corrida = dict(config)
        corrida["entorno"] = entorno
        observer = SeguimientoMetricasDecision(
            grilla, float(dossier["profundidad_s25_m"]), tuple(dossier["caras_s25"]))
        resultado, bordes, clima = demo.correr(
            corrida, grilla, vertidas, observadores=(observer,), modelo_hidratacion=fuente)
        if clima is None or clima.es_sintetica:
            raise ValueError("La corrida perdió el snapshot GEFS real")
        metricas = observer.resumen()
        fila = {
            "proteccion": etiqueta,
            "espesor_manta_m": (0.03 if "3cm" in etiqueta else
                0.05 if "5cm" in etiqueta else 0.06 if "6cm" in etiqueta else
                0.075 if "7_5cm" in etiqueta else 0.0),
            "espesor_aislacion_inferior_m": (0.05 if "base_5cm" in etiqueta else
                0.10 if "base_10cm" in etiqueta else 0.0),
            "retiro_manta_h_desde_colocacion": (120 if etiqueta == "manta_superior_3cm_120h" else
                168 if etiqueta == "manta_lateral_superior_3cm_168h" else
                360 if etiqueta.startswith("manta_5cm_360h_") else
                240 if etiqueta.startswith("manta_lateral_superior_5cm_240h") or
                      etiqueta.startswith("manta_5cm_240h_base_") else "actual"),
            "M_T_c": metricas["M_T_c"], "t_M_T_h": metricas["t_M_T_h"],
            "M_delta_c": metricas["M_delta_c"], "t_M_delta_h": metricas["t_M_delta_h"],
            "xyz_M_T_m": metricas["xyz_M_T_m"],
            "xyz_S25_min_m": metricas["xyz_S25_min_m"],
            "cumple_T_70c": metricas["M_T_c"] <= float(dossier["limites"]["temperatura_maxima_c"]),
            "cumple_delta_20c": metricas["M_delta_c"] <= float(dossier["limites"]["diferencia_interior_superficie_c"]),
            "residuo_balance_j": resultado.error_balance_j,
        }
        filas.append(fila)
        print(f"[{indice}/{len(opciones)}] {etiqueta}: M_T={fila['M_T_c']:.2f} °C, "
              f"M_delta={fila['M_delta_c']:.2f} °C, balance={fila['residuo_balance_j']:.2e} J")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as archivo:
        writer = csv.DictWriter(archivo, fieldnames=list(filas[0]))
        writer.writeheader()
        writer.writerows(filas)
    print(f"\nClima: {climate_cfg['id']} · material: {curva['id']} · plan: {plan.plan_id}")
    print(f"Tabla: {OUTPUT}")
    print("Moldaje lateral IFC: se mantiene a 72 h desde colocación en todas las opciones.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
