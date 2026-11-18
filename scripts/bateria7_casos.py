"""Batería de 7 casos × 8 curvas × 3 climas (336 h, suelo 2 m, malla 0,25 m), en paralelo.

    python scripts/bateria7_casos.py [--procs 8] [--solo-nuevos]

Usa exactamente la misma ruta de cálculo que ``evaluar_dossier_final.py`` (mismo solver,
entradas y observador S25); solo cambia que las corridas se reparten entre procesos y que
se guarda la serie temporal por paso. Configuración: ``configs/defensa/bateria7_casos.json``.

* Casos nuevos (2, 3, 6): 24 escenarios cada uno.
* Casos reutilizables (1, 4, 5, 7): solo 3 escenarios, para (a) comprobar que reproducen las
  métricas ya publicadas y (b) obtener series y campos 3D comparables.
* Campos 3D (float32, cada 6 h y en los instantes de los máximos) solo para 3 escenarios
  fijos a priori: M1_SS central, M4_SS frío y M4_SS central.

Resultados: ``results/bateria7/cache/<caso>__<escenario>.json`` (reanudable),
``results/bateria7/campos/*.npz``.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts"))
sys.path.insert(0, str(RAIZ / "src"))

CONFIG = RAIZ / "configs/defensa/bateria7_casos.json"
SAL = RAIZ / "results/bateria7"
NUEVOS = ("bat2_base_simultaneo", "bat3_base_largo_24h", "bat6_plan_largo_24h", "bat4_T11", "bat4_T12", "bat4_T13", "bat11_base_original_T10")
REUTIL = ("bat1_base_original", "bat4_plan1_original_10h", "bat5_plan_simultaneo", "bat7_control_sin_manta")
ESC_CAMPOS = ("M1_SS__GEFS_20261116_central", "M4_SS__GEFS_20261116_frio", "M4_SS__GEFS_20261116_central")
_G = {}


def _contexto():
    if _G:
        return _G
    import evaluar_dossier_final as ev
    import demo_pipeline_ifc_3d as pipe
    from mcplan.geometry.ifc import extraer_vertidas
    config, dossier = ev.cargar_configuracion_dossier(CONFIG)
    curvas = ev.leer_curvas(dossier)
    climas = dossier["climas"]
    matriz = ev.construir_matriz(curvas, climas, dossier)
    t0 = dossier["inicio_obra_utc"]
    config["tiempo_final_h"] = float(dossier.get("tiempo_final_h", 336.0))
    config["paso_voxel_m"] = float(dossier.get("paso_voxel_m", 0.25))
    config["paso_tiempo_max_s"] = float(dossier.get("paso_tiempo_max_s", 600.0))
    overlays = {str(c["id"]): ev.leer_jsonc((RAIZ / c["config_jsonc"]).resolve()) for c in climas}
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    planes_cfg = {str(p["id"]): p for p in dossier["planes"]}
    _G.update(ev=ev, config=config, dossier=dossier, matriz={m[0]: m for m in matriz}, t0=t0,
              overlays=overlays, originales=originales, planes_cfg=planes_cfg, grillas={})
    return _G


def correr_una(args):
    plan_id, esc_id = args
    destino = SAL / "cache" / f"{plan_id}__{esc_id}.json"
    if destino.exists():
        return json.loads(destino.read_text(encoding="utf-8"))
    g = _contexto()
    import demo_e10_entorno as demo
    from mcplan.geometry.voxelizacion import voxelizar_vertidas
    from mcplan.thermal.decision import SeguimientoMetricasDecision
    ev, dossier = g["ev"], g["dossier"]
    pcfg = g["planes_cfg"][plan_id]
    plan = ev.cargar_plan(pcfg, g["originales"])
    if plan_id not in g["grillas"]:
        vert = plan.aplicar(g["originales"])
        g["grillas"][plan_id] = (vert, voxelizar_vertidas(vert, g["config"]["paso_voxel_m"]))
    vertidas, grilla = g["grillas"][plan_id]
    _, (curva_cfg, hidr), clima_cfg, peso = g["matriz"][esc_id]
    cfg = dict(g["config"])
    ov = g["overlays"][str(clima_cfg["id"])]
    for campo in ("clima", "sitio", "entorno", "condicion_superficie"):
        if campo in ov:
            cfg[campo] = ov[campo]
    if pcfg.get("entorno") is not None:
        cfg["entorno"] = pcfg["entorno"]
    cfg.setdefault("clima", {})["inicio_obra_utc"] = g["t0"]
    obs = SeguimientoMetricasDecision(grilla, float(dossier.get("profundidad_s25_m", 0.025)),
                                      tuple(dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"])),
                                      guardar_historial=True)
    t_ini = time.time()
    res, _, aire = demo.correr(cfg, grilla, vertidas, observadores=(obs,), modelo_hidratacion=hidr)
    if aire is None or aire.es_sintetica:
        raise ValueError("clima no real")
    m = obs.resumen()
    # Serie por paso, reducida a ~10 min.
    h = obs.historial
    serie = {k: [round(x[k], 4) for x in h] for k in h[0]}
    out = {"plan_id": plan_id, "escenario_id": esc_id, "peso": peso, "metricas": m,
           "D_h": plan.duracion_obra_h(0.0), "segundos": time.time() - t_ini, "serie": serie}
    if esc_id in ESC_CAMPOS:
        vid = grilla.vertida_id
        solido = vid > 0
        horas = np.asarray(res.tiempos_s) / 3600.0
        idx = sorted({int(round(x)) for x in (*np.arange(0, 337, 6.0), m["t_M_T_h"], m["t_M_delta_h"])
                      if 0 <= round(x) < len(horas)})
        campos = np.stack([np.asarray(res.temperaturas_c[i])[solido] for i in idx]).astype(np.float32)
        (SAL / "campos").mkdir(parents=True, exist_ok=True)
        np.savez_compressed(SAL / "campos" / f"{plan_id}__{esc_id}.npz", T=campos, t_h=horas[idx].astype(np.float32),
                            vid=vid[solido].astype(np.uint8), forma=np.asarray(grilla.forma),
                            origen=np.asarray(grilla.origen_m), paso=np.asarray(grilla.paso_m),
                            mascara=solido, ids_vertida=vid.astype(np.uint8))
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(json.dumps(out), encoding="utf-8")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--solo-nuevos", action="store_true")
    a = ap.parse_args()
    g = _contexto()
    esc = list(g["matriz"])
    tareas = [(p, e) for p in NUEVOS for e in esc]
    if not a.solo_nuevos:
        tareas += [(p, e) for p in REUTIL for e in ESC_CAMPOS]
    print(f"{len(tareas)} corridas, {a.procs} procesos", flush=True)
    t0 = time.time()
    ok = 0
    with mp.get_context("spawn").Pool(a.procs) as pool:
        for r in pool.imap_unordered(correr_una, tareas):
            ok += 1
            print(f"[{ok}/{len(tareas)}] {r['plan_id']} {r['escenario_id']} M_T={r['metricas']['M_T_c']:.2f} "
                  f"M_d={r['metricas']['M_delta_c']:.3f} ({r['segundos']:.0f} s, total {time.time() - t0:.0f} s)", flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
