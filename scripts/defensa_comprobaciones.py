"""Comprobaciones posteriores al informe: colocación progresiva a caudal constante,
reglas de retiro escalonado de la manta y variantes sin aislación inferior.

    python scripts/defensa_comprobaciones.py [--procs 15] [grupo ...]    # grupos: prog retiro sinAI
Resultados en results/defensa/cache2/*.json (reanudable) y tabla results/defensa/comprobaciones.csv.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "scripts"))
CACHE2 = RAIZ / "results/defensa/cache2"
CAPAS_N = 4   # la manta de 80 mm se trata como 4 capas de 20 mm retirables por separado


def planes_base() -> dict:
    import defensa_busqueda as db
    out = {p["id"]: p for p in db._leer_jsonc(RAIZ / "configs/defensa/planes_final.json")["planes"]}
    return out


def retiro_capas(spec, inicio_h, intervalo_h, n=CAPAS_N):
    """Retiro escalonado: capa j (j=0..n-1) se retira en inicio_h + j·intervalo_h desde t0."""
    s = copy.deepcopy(spec)
    esp = s["manta"]["espesor_m"]
    s["manta_capas"] = {"n": n, "espesor_capa_m": esp / n, "k": s["manta"].get("k", 0.025),
                        "inicio_h": s["manta"].get("inicio_h", 15.0),
                        "retiros_h": None if inicio_h is None else [inicio_h + j * intervalo_h for j in range(n)]}
    s.pop("manta")
    return s


def retiro_lista(spec, dias, n=CAPAS_N):
    """Retiro explícito: la capa j se retira el día dias[j] desde t0 (manta de 80 mm = 4 capas de 20 mm)."""
    s = retiro_capas(spec, 0.0, 0.0, n)
    s["manta_capas"]["retiros_h"] = [24.0 * d for d in dias]
    return s


REGLAS_RET2 = {   # días desde t0 de la retirada de cada capa de 20 mm
    "A_14d_desde56_(falla)": (56, 70, 84, 98),
    "B_21d_desde56": (56, 77, 98, 119),
    "C_28d_desde56": (56, 84, 112, 140),
    "D_14d_ultima_d140": (56, 70, 84, 140),
    "E_14d_ultima_d168": (56, 70, 84, 168),
    "F_21d_desde84": (84, 105, 126, 147),
    "G_28d_ultima_d168": (56, 84, 112, 168),
}


def tareas() -> list[dict]:
    import defensa_generar_planes as gen
    base = planes_base()
    U63, U36 = base["f_P2_L2m_U63"], base["f_P2_L2m_U36"]
    esc = (("M4_SS", "fria"), ("M4_SS", "calida"), ("M1_SS", "central"))
    t = []
    for spec in (U63, U36):                                   # 2. colocación progresiva
        for c, w in esc:
            t.append({"grupo": "prog", "nombre": f"prog_{spec['id']}", "spec": spec, "curva": c, "clima": w,
                      "horizonte_h": 1344.0, "progresiva": True})
    reglas = {"base_sin_retiro": None, "R56_c14d": (1344.0, 336.0), "R42_c7d": (1008.0, 168.0),
              "R56_c7d": (1344.0, 168.0), "R84_c14d": (2016.0, 336.0)}
    for nombre, regla in reglas.items():                       # 3. reglas de retiro
        spec = retiro_capas(U63, *(regla if regla else (None, 0)))
        spec["id"] = f"ret_{nombre}"
        for c, w in esc:
            t.append({"grupo": "retiro", "nombre": spec["id"], "spec": spec, "curva": c, "clima": w,
                      "horizonte_h": 3600.0, "progresiva": False})
    for tc, mm, suelo in [(10.0, m, s) for m in (100, 120, 150) for s in ("L2m", "L0_8m")] + \
                         [(15.0, m, s) for m in (120, 150) for s in ("L2m", "L0_8m")]:   # 4. sin aislación inferior
        pid = f"sAI_T{tc:.0f}_PUR{mm}_{suelo}"
        spec = gen.spec_plan(pid, "sin_aislacion_inferior", gen.paquete("x", tcol=tc, manta_mm=mm, aisl_mm=0, suelo=suelo),
                             tasa=gen.R_ALTA)
        t.append({"grupo": "sinAI", "nombre": pid, "spec": spec, "curva": "M4_SS", "clima": "fria",
                  "horizonte_h": 1344.0, "progresiva": False})
    U36b = copy.deepcopy(U36); U36b["id"] = "f_P2_L0_8m_U36"; U36b["suelo"] = dict(U36["suelo"], L=0.8)
    for spec in (U36, U36b):                                  # candidato: evaluación completa 24 escenarios
        for c in ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"):
            for w in ("fria", "central", "calida"):
                t.append({"grupo": "cand", "nombre": f"cand_{spec['id']}", "spec": spec, "curva": c, "clima": w,
                          "horizonte_h": 1344.0, "progresiva": True})
    for spec in (U36, U36b):                                  # retiro 20 mm c/14 d desde d56, candidato completo
        sr = retiro_capas(spec, 1344.0, 336.0)
        sr["id"] = f"retC_{spec['id']}"
        for c in ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"):
            for w in ("fria", "central", "calida"):
                t.append({"grupo": "retC", "nombre": sr["id"], "spec": sr, "curva": c, "clima": w,
                          "horizonte_h": 3600.0, "progresiva": True})
    # --- Candidato con aislación inferior de 100 mm (muestreo reducido) ---
    def con_ai(spec, nombre, mm=100):
        s2 = copy.deepcopy(spec); s2["id"] = nombre
        s2["aislacion_inferior"] = dict(s2["aislacion_inferior"], espesor_m=mm / 1000.0)
        return s2
    c24 = [(c, w) for c in ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A")
           for w in ("fria", "central", "calida")]
    for base_spec, suf in ((U36, "L2m"), (U36b, "L0_8m")):
        a = con_ai(base_spec, f"ai100_{suf}")
        r = retiro_capas(a, 1344.0, 336.0); r["id"] = f"ai100ret_{suf}"
        for c, w in c24:
            t.append({"grupo": "ai100", "nombre": a["id"], "spec": a, "curva": c, "clima": w, "horizonte_h": 1344.0,
                      "progresiva": True, "reducido": True})
            t.append({"grupo": "ai100ret", "nombre": r["id"], "spec": r, "curva": c, "clima": w, "horizonte_h": 3600.0,
                      "progresiva": True, "reducido": True})
    for base_spec, suf in ((U36, "L2m"), (U36b, "L0_8m")):    # reglas de retiro más gradual/tardío, M4_SS cálida
        a = con_ai(base_spec, f"ai100_{suf}")
        for nom, dias in REGLAS_RET2.items():
            r = retiro_lista(a, dias); r["id"] = f"ret2_{nom}_{suf}"
            t.append({"grupo": "ret2", "nombre": r["id"], "spec": r, "curva": "M4_SS", "clima": "calida",
                      "horizonte_h": 5400.0, "progresiva": True, "reducido": True})
    for base_spec, suf in ((U36, "L2m"), (U36b, "L0_8m")):    # regla B, matriz completa
        a = con_ai(base_spec, f"ai100_{suf}")
        r = retiro_lista(a, (56, 77, 98, 119)); r["id"] = f"retB_{suf}"
        for c, w in c24:
            t.append({"grupo": "retB", "nombre": r["id"], "spec": r, "curva": c, "clima": w,
                      "horizonte_h": 3400.0, "progresiva": True, "reducido": True})
    a0 = con_ai(U36b, "ai100dt_L0_8m")
    for dt in (600.0, 300.0, 150.0):
        t.append({"grupo": "dt", "nombre": f"ai100dt{int(dt)}_L0_8m", "spec": a0, "curva": "M4_SS", "clima": "calida",
                  "horizonte_h": 1344.0, "progresiva": True, "reducido": True, "dt_max_s": dt})
    for spec, esc2 in ((U36, esc), (U36b, esc[:1])):          # duraciones según el volumen numérico de cada vertida
        sv = copy.deepcopy(spec); sv["vol_numerico"] = True
        for c, w in esc2:
            t.append({"grupo": "vol", "nombre": f"vol_{spec['id']}", "spec": sv, "curva": c, "clima": w,
                      "horizonte_h": 1344.0, "progresiva": True})
    for x in t:
        clave = json.dumps({k: x[k] for k in ("nombre", "spec", "curva", "clima", "horizonte_h", "progresiva", "reducido", "dt_max_s") if k in x}, sort_keys=True)
        x["archivo"] = str(CACHE2 / f"{x['nombre']}__{x['curva']}__{x['clima']}__{hashlib.sha256(clave.encode()).hexdigest()[:10]}.json")
    return t


def tiempos_frente(grilla, vertidas, hor, paso_s=300.0):
    """Colocación progresiva a caudal constante: cada vertida se llena por capas (z ascendente) y,
    dentro de la capa, barriendo en x e y. Los voxeles tienen igual volumen, así que el instante
    del voxel k de N es inicio + duración·k/N (caudal = volumen/duración, constante). Cuantizado a 5 min."""
    tiempos = np.full(grilla.forma, np.inf)
    for v in vertidas:
        i = int(v.vertida_id)
        ijk = np.argwhere(grilla.vertida_id == i)
        orden = np.lexsort((ijk[:, 1], ijk[:, 0], ijk[:, 2]))
        n = len(ijk)
        frac = (np.arange(n) + 1.0) / n
        t = (hor["activaciones_h"][i] + hor["duraciones_h"][i] * frac) * 3600.0
        t = np.round(t / paso_s) * paso_s
        t = np.minimum(t, (hor["activaciones_h"][i] + hor["duraciones_h"][i]) * 3600.0)
        sel = ijk[orden]
        tiempos[sel[:, 0], sel[:, 1], sel[:, 2]] = t
        # el primer instante activo no puede ser anterior al inicio de la vertida
        tiempos[sel[:, 0], sel[:, 1], sel[:, 2]] = np.maximum(tiempos[sel[:, 0], sel[:, 1], sel[:, 2]],
                                                              hor["activaciones_h"][i] * 3600.0)
    return tiempos


def ejecutar(x: dict) -> dict:
    import defensa_busqueda as db
    import demo_e10_entorno as demo
    from mcplan.geometry.ifc import extraer_vertidas
    from mcplan.geometry.voxelizacion import voxelizar_vertidas
    from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision
    from mcplan.thermal.hidratacion import ModeloFuenteTabulada
    t_ini = time.time()
    spec = copy.deepcopy(x["spec"])
    if spec.get("vol_numerico"):
        o0 = extraer_vertidas((RAIZ / "models/ifc/Fundacion_1_con_propiedades.ifc").resolve())
        g0 = voxelizar_vertidas(o0, 0.25)
        vol = {i: float((g0.vertida_id == i).sum() * np.prod(g0.paso_m)) for i in range(1, 6)}
        tasa, tt, act, dur = float(spec["tasa_m3_h"]), 0.0, {}, {}
        for i in range(1, 6):
            act[str(i)], dur[str(i)] = tt, vol[i] / tasa
            tt += dur[str(i)]
        spec.pop("tasa_m3_h"); spec.pop("espera_h", None); spec.pop("inicio_primera_h", None)
        spec["activaciones_h"], spec["duraciones_h"] = act, dur
    hor = db.construir_plan(spec)
    capas = spec.pop("manta_capas", None)
    entorno = db.construir_entorno(spec, hor)
    if capas:
        for i in range(1, 6):
            lista = entorno["protecciones_por_vertida"].setdefault(str(i), [])
            for j in range(capas["n"]):
                ret = None if capas["retiros_h"] is None else capas["retiros_h"][j]
                lista.append({"nombre": f"manta_capa{j + 1}", "caras": ["x-", "x+", "y-", "y+", "z+"],
                              "espesor_m": capas["espesor_capa_m"], "conductividad_w_m_k": capas["k"],
                              "inicio_h": max(0.0, capas["inicio_h"] - hor["activaciones_h"][i]),
                              "inicio_desde": "colocacion", "retiro_h": ret, "retiro_desde": "inicio_obra"})
    config = {}
    for capa in ("configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"):
        config.update(db._leer_jsonc(RAIZ / capa))
    config.pop("sobrescribir_vertidas", None)
    config.update({"tiempo_final_h": x["horizonte_h"], "paso_voxel_m": 0.25,
                   "paso_tiempo_max_s": float(x.get("dt_max_s", 600.0)), "entorno": entorno})
    config["clima"] = dict(db._leer_jsonc(RAIZ / db.CLIMAS[x["clima"]][0])["clima"])
    config["clima"]["inicio_obra_utc"] = "2026-11-19T13:00:00Z"
    if x["horizonte_h"] > 672.0:
        config["clima"]["snapshot"] = dict(config["clima"]["snapshot"], extension="ciclo_diario")
    config["sitio"] = {"grupo": 3, "latitud": -31.1115, "longitud": -71.6147}
    tcol = spec.get("temperatura_colocacion_c", 18.0)
    originales = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan = PlanConstruccion(spec["id"], hor["activaciones_h"], hor["duraciones_h"], {i: float(tcol) for i in range(1, 6)},
                            frentes_simultaneos_max=int(spec.get("frentes_max", 1)))
    vertidas = plan.aplicar(originales)
    grilla = voxelizar_vertidas(vertidas, 0.25)
    tiempos = tiempos_frente(grilla, vertidas, hor) if x["progresiva"] else None
    caras = ("x-", "x+", "y-", "y+", "z-", "z+")
    o_fin = SeguimientoMetricasDecision(grilla, 0.025, caras)
    o14 = SeguimientoMetricasDecision(grilla, 0.025, caras, tiempo_max_s=336 * 3600.0)
    o56 = SeguimientoMetricasDecision(grilla, 0.025, caras, tiempo_max_s=1344 * 3600.0)

    class Serie:
        def __init__(self):
            self.h, self.d, self.tm, self.s = [], [], [], []
            self.sig = 0.0

        def registrar(self, *a):
            pass

        def registrar_paso(self, t_s, dt_s, solver, T, f):
            u = o_fin.ultimo
            if u and t_s + 1e-6 >= self.sig and abs(u["t_s"] - t_s) < 1e-6:
                self.h.append(round(t_s / 3600, 3)); self.d.append(u["delta_c"]); self.tm.append(u["T_max_c"])
                self.s.append(u["S25_min_c"]); self.sig = t_s + 3 * 3600.0

    serie = Serie()
    fuente = ModeloFuenteTabulada.desde_source_terms(db.SOURCE_TERMS, x["curva"], cola_cero=True)
    if x.get("reducido"):
        import defensa_malla_candidato as dm
        res, _, _ = dm.correr_reducido(config, grilla, vertidas, (o_fin, o14, o56, serie), fuente, tiempos, ())
    else:
        res, _, _ = demo.correr(config, grilla, vertidas, observadores=(o_fin, o14, o56, serie),
                                modelo_hidratacion=fuente, tiempos_activacion_s_por_voxel=tiempos)
    out = {"nombre": x["nombre"], "grupo": x["grupo"], "curva": x["curva"], "clima": x["clima"],
           "progresiva": x["progresiva"], "horizonte_h": x["horizonte_h"], "D_h": hor["D_h"],
           "retiros_h": (spec.get("manta_capas") or capas or {}).get("retiros_h") if capas else None,
           "h14": o14.resumen(), "h56": o56.resumen(), "hfin": o_fin.resumen(),
           "serie": {"h": serie.h, "delta_c": serie.d, "Tmax_c": serie.tm, "S25_min_c": serie.s},
           "balance_error_J": float(res.error_balance_j), "segundos": time.time() - t_ini}
    Path(x["archivo"]).write_text(json.dumps(out), encoding="utf-8")
    return {"n": x["nombre"], "c": x["curva"], "w": x["clima"], "d56": out["h56"]["M_delta_c"],
            "dfin": out["hfin"]["M_delta_c"], "tfin": out["hfin"]["t_M_delta_h"], "s": out["segundos"]}


def _w(x):
    try:
        return ejecutar(x)
    except Exception as e:  # noqa: BLE001
        return {"error": f"{type(e).__name__}: {e}", "n": x["nombre"], "c": x["curva"], "w": x["clima"]}


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = sys.argv[1:]
    procs = 15
    if "--procs" in args:
        k = args.index("--procs"); procs = int(args[k + 1]); del args[k:k + 2]
    grupos = set(args) or {"prog", "retiro", "sinAI"}
    CACHE2.mkdir(parents=True, exist_ok=True)
    ts = [x for x in tareas() if x["grupo"] in grupos]
    pend = [x for x in ts if not Path(x["archivo"]).exists()]
    print(f"{len(ts)} corridas, {len(pend)} por ejecutar", flush=True)
    if pend:
        with mp.get_context("spawn").Pool(min(procs, len(pend))) as pool:
            for k, r in enumerate(pool.imap_unordered(_w, pend), 1):
                print(f"[{k}/{len(pend)}] {r}", flush=True)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
