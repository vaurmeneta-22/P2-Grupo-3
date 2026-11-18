"""Búsqueda de planes para la defensa: corre plan × curva × clima con caché reanudable.

Cada corrida (una combinación plan/curva/clima) se guarda como un JSON en
``results/defensa/cache`` con sus métricas en los horizontes de 14 y 28 días, la
serie horaria Tmax/S25/ΔT y el balance de energía. Si el archivo existe y las entradas
(plan, curva, clima, malla, paso, horizonte y código) no cambiaron, no se repite.

Uso (desde la raíz, con ``PYTHONPATH=src``)::

    python scripts/defensa_busqueda.py configs/defensa/planes_cribado.json \
        --curvas M1_SS --climas fria,central,calida --procs 8

El plan se declara con una tasa de colocación y esperas (ver ``construir_plan``);
la activación térmica de cada vertida es instantánea al inicio de su colocación y D
usa el FIN físico de colocación.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "scripts"))

CACHE = RAIZ / "results/defensa/cache"
SOURCE_TERMS = RAIZ / "data/hidratacion/source_terms.txt"
CURVAS = ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A")
CLIMAS = {  # pesos = fracción de los 31 miembros GEFS que representa cada trayectoria
    "fria": ("configs/climas/gefs_20261116_28d/fria.jsonc", 11 / 31),
    "central": ("configs/climas/gefs_20261116_28d/central.jsonc", 10 / 31),
    "calida": ("configs/climas/gefs_20261116_28d/calida.jsonc", 10 / 31),
}
# Volúmenes por vertida del IFC (docs/dossier_final.md), m³.
VOLUMENES_M3 = {1: 147.626, 2: 147.626, 3: 214.224, 4: 88.925, 5: 38.375}
ARCHIVOS_CODIGO = ("src/mcplan/thermal/decision.py", "src/mcplan/thermal/hidratacion.py",
                   "src/mcplan/thermal/vertidas.py", "src/mcplan/thermal/bordes.py",
                   "src/mcplan/thermal/diffusion_nd.py", "scripts/demo_e10_entorno.py",
                   "scripts/defensa_busqueda.py")
_GRILLAS: dict = {}


def _leer_jsonc(ruta: Path) -> dict:
    import demo_pipeline_ifc_3d as pipe
    return json.loads(pipe.quitar_comentarios_jsonc(ruta.read_text(encoding="utf-8")))


def _sha(ruta: Path) -> str:
    return hashlib.sha256(ruta.read_bytes()).hexdigest()


def huella_codigo() -> str:
    h = hashlib.sha256()
    for r in ARCHIVOS_CODIGO:
        h.update(_sha(RAIZ / r).encode())
    return h.hexdigest()[:12]


# ----------------------------------------------------------------------------------
def construir_plan(spec: dict) -> dict:
    """Resuelve horarios: devuelve activaciones, duraciones, D y tasa efectiva.

    Campos: ``tasa_m3_h`` (una sola cuadrilla, sin traslapes), ``espera_h`` (escalar o lista
    de 4 esperas entre vertidas), ``inicio_primera_h`` (postergación de la primera colocación
    desde t0), o bien ``activaciones_h`` y ``duraciones_h`` explícitas.
    """
    if "activaciones_h" in spec:
        act = {int(k): float(v) for k, v in spec["activaciones_h"].items()}
        dur = {int(k): float(v) for k, v in spec["duraciones_h"].items()}
    else:
        tasa = float(spec["tasa_m3_h"])
        esperas = spec.get("espera_h", 0.0)
        esperas = list(esperas) if isinstance(esperas, (list, tuple)) else [float(esperas)] * 4
        t = float(spec.get("inicio_primera_h", 0.0))
        act, dur = {}, {}
        for i in range(1, 6):
            act[i], dur[i] = t, VOLUMENES_M3[i] / tasa
            t += dur[i] + (esperas[i - 1] if i < 5 else 0.0)
    fin = max(act[i] + dur[i] for i in act)
    solape = any(act[i + 1] < act[i] + dur[i] - 1e-9 for i in range(1, 5))
    if solape and int(spec.get("frentes_max", 1)) < 2:
        raise ValueError(f"{spec['id']}: vertidas solapadas con un solo frente")
    return {"activaciones_h": act, "duraciones_h": dur, "D_h": fin, "t0_h": 0.0,
            "tasa_media_m3_h": sum(VOLUMENES_M3.values()) / max(fin - min(act.values()), 1e-9)}


def construir_entorno(spec: dict, horarios: dict) -> dict:
    mold = spec.get("moldaje", {"vertidas": [1, 2, 5], "espesor_m": 0.018, "k": 0.1})
    manta = spec.get("manta")
    aisl = spec.get("aislacion_inferior")
    suelo = spec.get("suelo", {"T": 20.0, "k": 1.5, "L": 2.0})
    prot = {}
    for i in range(1, 6):
        lista = []
        if mold and i in mold.get("vertidas", []):
            lista.append({"nombre": "moldaje", "caras": ["x-", "x+", "y-", "y+"],
                          "espesor_m": mold["espesor_m"], "conductividad_w_m_k": mold["k"],
                          "retiro_h": mold.get("retiro_h", "ifc"), "retiro_desde": "colocacion"})
        if manta:
            inicio_rel = max(0.0, float(manta.get("inicio_h", 15.0)) - horarios["activaciones_h"][i])
            lista.append({"nombre": "manta", "caras": ["x-", "x+", "y-", "y+", "z+"],
                          "espesor_m": manta["espesor_m"], "conductividad_w_m_k": manta.get("k", 0.025),
                          "inicio_h": inicio_rel, "inicio_desde": "colocacion",
                          "retiro_h": manta.get("retiro_h"), "retiro_desde": "instalacion"})
        if aisl and i == 1:
            lista.append({"nombre": "aislacion_inferior", "caras": ["z-"], "espesor_m": aisl["espesor_m"],
                          "conductividad_w_m_k": aisl.get("k", 0.04), "retiro_h": None,
                          "retiro_desde": "colocacion"})
        if lista:
            prot[str(i)] = lista
    return {"activo": True, "aire": {"h_w_m2_k": float(spec.get("h_aire", 8.0))},
            "suelo": {"activo": True, "temperatura_c": suelo["T"], "conductividad_w_m_k": suelo["k"],
                      "espesor_efectivo_m": suelo["L"]},
            "protecciones_por_vertida": prot}


def tarea(spec: dict, curva: str, clima: str, horizonte_h: float, defecto: dict) -> dict:
    t = {"plan": spec, "curva": curva, "clima": clima, "horizonte_h": horizonte_h,
         "dx_m": float(spec.get("dx_m", defecto["dx_m"])),
         "dt_max_s": float(spec.get("dt_max_s", defecto["dt_max_s"])),
         "s25_m": float(spec.get("s25_m", defecto["s25_m"])),
         "guardar_campos_h": spec.get("guardar_campos_h", [])}
    clave = json.dumps(t, sort_keys=True) + huella_codigo() + _sha(SOURCE_TERMS)[:12] + \
        _sha(RAIZ / CLIMAS[clima][0].replace(".jsonc", ".jsonc"))[:12]
    t["clave"] = hashlib.sha256(clave.encode()).hexdigest()[:16]
    t["archivo"] = str(CACHE / f"{spec['id']}__{curva}__{clima}__{t['clave']}.json")
    return t


# ----------------------------------------------------------------------------------
def ejecutar(t: dict) -> dict:
    """Corre una combinación. Se ejecuta en un proceso hijo."""
    import demo_e10_entorno as demo
    from mcplan.geometry.ifc import extraer_vertidas
    from mcplan.geometry.voxelizacion import voxelizar_vertidas
    from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision
    from mcplan.thermal.hidratacion import ModeloFuenteTabulada

    t_ini = time.time()
    spec = t["plan"]
    hor = construir_plan(spec)
    config = {}
    for capa in ("configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"):
        config.update(_leer_jsonc(RAIZ / capa))
    config.pop("sobrescribir_vertidas", None)
    config.update({"tiempo_final_h": t["horizonte_h"], "paso_voxel_m": t["dx_m"],
                   "paso_tiempo_max_s": t["dt_max_s"], "entorno": construir_entorno(spec, hor)})
    overlay = _leer_jsonc(RAIZ / CLIMAS[t["clima"]][0])
    config["clima"] = dict(overlay["clima"])
    config["clima"]["inicio_obra_utc"] = "2026-11-19T13:00:00Z"
    if t["horizonte_h"] > 672.0:  # continuación etiquetada: repite el último ciclo diario del pronóstico
        config["clima"]["snapshot"] = dict(config["clima"]["snapshot"], extension="ciclo_diario")
    config["sitio"] = {"grupo": 3, "latitud": -31.1115, "longitud": -71.6147}
    tcol = spec.get("temperatura_colocacion_c", 18.0)
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan = PlanConstruccion(spec["id"], hor["activaciones_h"], hor["duraciones_h"],
                            {i: float(tcol[str(i)] if isinstance(tcol, dict) else tcol) for i in range(1, 6)},
                            frentes_simultaneos_max=int(spec.get("frentes_max", 1)))
    vertidas = plan.aplicar(originales)
    dx = t["dx_m"]
    if dx not in _GRILLAS:
        _GRILLAS[dx] = voxelizar_vertidas(vertidas, dx)
    grilla = _GRILLAS[dx]
    fuente = ModeloFuenteTabulada.desde_source_terms(SOURCE_TERMS, t["curva"], cola_cero=True)
    caras = ("x-", "x+", "y-", "y+", "z-", "z+")
    obs_total = SeguimientoMetricasDecision(grilla, t["s25_m"], caras)
    obs_14 = SeguimientoMetricasDecision(grilla, t["s25_m"], caras, tiempo_max_s=336.0 * 3600.0)
    obs_28 = SeguimientoMetricasDecision(grilla, t["s25_m"], caras, tiempo_max_s=672.0 * 3600.0)
    # Sensibilidad de definición: S25 sin la cara inferior (contacto con el suelo). No sustituye a M_delta.
    obs_sin_base = SeguimientoMetricasDecision(grilla, t["s25_m"], ("x-", "x+", "y-", "y+", "z+"))

    class Serie:
        def __init__(self):
            self.h, self.tmax, self.s25, self.delta = [], [], [], []
            self._sig = 0.0

        def registrar(self, tiempo_s, temperatura_c, solver):
            pass  # el solver exige esta interfaz; la serie se llena en registrar_paso

        def registrar_paso(self, tiempo_s, dt_s, solver, temperatura_c, fuente_w_m3):
            u = obs_total.ultimo
            if u and tiempo_s + 1e-6 >= self._sig and abs(u["t_s"] - tiempo_s) < 1e-6:
                self.h.append(round(tiempo_s / 3600.0, 4)); self.tmax.append(u["T_max_c"])
                self.s25.append(u["S25_min_c"]); self.delta.append(u["delta_c"])
                self._sig = tiempo_s + 3600.0

    serie = Serie()
    guardar = [float(x) for x in t["guardar_campos_h"]]
    resultado, bordes, aire = demo.correr(config, grilla, vertidas,
                                          observadores=(obs_total, obs_14, obs_28, obs_sin_base, serie),
                                          modelo_hidratacion=fuente)
    if aire is None or aire.es_sintetica:
        raise ValueError("clima sintético o ausente")
    mfin, m28, m14 = obs_total.resumen(), obs_28.resumen(), obs_14.resumen()
    out = {k: t[k] for k in ("curva", "clima", "horizonte_h", "dx_m", "dt_max_s", "s25_m", "clave")}
    out.update({
        "plan_id": spec["id"], "categoria": spec.get("categoria", ""), "plan": spec, "horarios": hor,
        "D_h": hor["D_h"], "peso_clima": CLIMAS[t["clima"]][1],
        "h14": m14, "h28": m28, "hfin": mfin, "hfin_sin_base": obs_sin_base.resumen(),
        "balance": {"error_J": float(resultado.error_balance_j), "energia_fuente_J": float(resultado.energia_fuente_j),
                    "energia_activacion_J": float(getattr(resultado, "energia_activacion_j", 0.0)),
                    "energia_frontera_J": float(sum(resultado.energia_frontera_j.values()))
                    if isinstance(getattr(resultado, "energia_frontera_j", None), dict) else None},
        "serie": {"h": serie.h, "Tmax_c": serie.tmax, "S25_min_c": serie.s25, "delta_c": serie.delta},
        "aire": {"fuente": aire.descripcion.get("fuente"), "snapshot": aire.descripcion.get("snapshot"),
                 "cobertura": aire.descripcion.get("cobertura")},
        "segundos": time.time() - t_ini, "codigo": huella_codigo(),
    })
    if guardar:
        tiempos = np.asarray(resultado.tiempos_s) / 3600.0
        campos = {}
        for h in guardar:
            j = int(np.argmin(np.abs(tiempos - h)))
            campos[f"T_{tiempos[j]:.0f}h"] = np.asarray(resultado.temperaturas_c[j], dtype=np.float32)
            campos[f"activa_{tiempos[j]:.0f}h"] = np.asarray(resultado.mascaras_activas[j], dtype=bool)
        np.savez_compressed(str(Path(t["archivo"]).with_suffix(".npz")), **campos,
                            vertida_id=np.asarray(grilla.vertida_id), origen=np.asarray(grilla.origen_m),
                            paso=np.asarray(grilla.paso_m))
    Path(t["archivo"]).write_text(json.dumps(out), encoding="utf-8")
    return {"archivo": t["archivo"], "plan": spec["id"], "curva": t["curva"], "clima": t["clima"],
            "M_T14": m14["M_T_c"], "M_d14": m14["M_delta_c"], "M_T28": m28["M_T_c"],
            "M_d28": m28["M_delta_c"], "seg": out["segundos"]}


def _worker(t):
    try:
        return ejecutar(t)
    except Exception as exc:  # noqa: BLE001 - se informa y se sigue con el resto
        return {"error": f"{type(exc).__name__}: {exc}", "plan": t["plan"]["id"], "curva": t["curva"], "clima": t["clima"]}


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("planes")
    ap.add_argument("--curvas", default="M1_SS", help="lista separada por comas o 'todas'")
    ap.add_argument("--climas", default="fria,central,calida")
    ap.add_argument("--horizonte-h", type=float, default=672.0)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--solo", default="", help="ids de plan separados por comas")
    ap.add_argument("--dx", type=float, default=0.25)
    ap.add_argument("--dt", type=float, default=600.0)
    ap.add_argument("--s25", type=float, default=0.025)
    a = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    planes = _leer_jsonc(Path(a.planes) if Path(a.planes).is_absolute() else RAIZ / a.planes)["planes"]
    if a.solo:
        planes = [p for p in planes if p["id"] in a.solo.split(",")]
    curvas = CURVAS if a.curvas == "todas" else tuple(a.curvas.split(","))
    climas = tuple(CLIMAS) if a.climas == "todas" else tuple(a.climas.split(","))
    defecto = {"dx_m": a.dx, "dt_max_s": a.dt, "s25_m": a.s25}
    tareas = [tarea(p, c, w, a.horizonte_h, defecto) for p in planes for c in curvas for w in climas]
    pendientes = [t for t in tareas if not Path(t["archivo"]).exists()]
    print(f"{len(tareas)} corridas, {len(tareas) - len(pendientes)} en caché, {len(pendientes)} por ejecutar; "
          f"procs={a.procs}; código={huella_codigo()}", flush=True)
    t0 = time.time()
    errores = 0
    if pendientes:
        with mp.get_context("spawn").Pool(min(a.procs, len(pendientes))) as pool:
            for k, r in enumerate(pool.imap_unordered(_worker, pendientes), 1):
                if "error" in r:
                    errores += 1
                    print(f"[{k}/{len(pendientes)}] ERROR {r['plan']} {r['curva']} {r['clima']}: {r['error']}", flush=True)
                else:
                    print(f"[{k}/{len(pendientes)}] {r['plan']} {r['curva']} {r['clima']}: "
                          f"M_T={r['M_T14']:.2f}/{r['M_T28']:.2f} M_Δ={r['M_d14']:.2f}/{r['M_d28']:.2f} "
                          f"(14d/28d) {r['seg']:.0f}s", flush=True)
    print(f"Listo en {time.time() - t0:.0f} s; errores={errores}")
    return 1 if errores else 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
