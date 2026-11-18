"""Comprobación espacial del candidato (36 m³/h progresivo, 10 °C, PUR 80 mm, AI 50 mm, suelo 0,8 m,
M4_SS cálida). Mallas: 'ref' (0,25³), 'iso125' (0,125³) y 'aniso005' (0,25×0,25×0,05).

    python scripts/defensa_malla_candidato.py verificar aniso005     # solo comprobaciones previas
    python scripts/defensa_malla_candidato.py correr iso125
    python scripts/defensa_malla_candidato.py correr aniso005
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import time

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "scripts"))

import defensa_busqueda as db  # noqa: E402
import defensa_comprobaciones as dc  # noqa: E402
import demo_e10_entorno as demo  # noqa: E402
import demo_pipeline_ifc_3d as pipe  # noqa: E402
from mcplan.geometry.ifc import extraer_vertidas  # noqa: E402
from mcplan.geometry.voxelizacion import voxelizar_vertidas  # noqa: E402
from mcplan.thermal.decision import PlanConstruccion, SeguimientoMetricasDecision  # noqa: E402
from mcplan.thermal.diffusion_nd import DifusionCartesiana  # noqa: E402
from mcplan.thermal.hidratacion import ModeloFuenteTabulada  # noqa: E402

MALLAS = {"ref": (0.25, 0.25, 0.25), "iso125": (0.125, 0.125, 0.125), "aniso005": (0.25, 0.25, 0.05),
          "comb": (0.125, 0.125, 0.05)}
CURVA, CLIMA, PLAN = "M4_SS", "calida", "f_P2_L0_8m_U36"
SALIDA = RAIZ / "results/defensa/cache2"
IFC = {1: 147.626, 2: 147.626, 3: 214.224, 4: 88.925, 5: 38.375}


REGLAS = {"B": (56, 77, 98, 119), "F": (84, 105, 126, 147), "A": (56, 70, 84, 98)}


def preparar(malla, ai_mm=50, suelo_L=0.8, regla=None, horizonte_h=1344.0):
    base = dc.planes_base()["f_P2_L2m_U36"]
    spec = copy.deepcopy(base); spec["id"] = PLAN; spec["suelo"] = dict(spec["suelo"], L=suelo_L)
    spec["aislacion_inferior"] = dict(spec["aislacion_inferior"], espesor_m=ai_mm / 1000.0)
    capas = None
    if regla:
        spec = dc.retiro_lista(spec, REGLAS[regla]); spec["id"] = PLAN
    hor = db.construir_plan(spec)
    capas = spec.pop("manta_capas", None)
    config = {}
    for capa in ("configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"):
        config.update(db._leer_jsonc(RAIZ / capa))
    config.pop("sobrescribir_vertidas", None)
    config.update({"tiempo_final_h": horizonte_h, "paso_voxel_m": MALLAS[malla][0], "paso_tiempo_max_s": 600.0,
                   "entorno": db.construir_entorno(spec, hor)})
    if capas:
        for i in range(1, 6):
            lista = config["entorno"]["protecciones_por_vertida"].setdefault(str(i), [])
            for j in range(capas["n"]):
                lista.append({"nombre": f"manta_capa{j + 1}", "caras": ["x-", "x+", "y-", "y+", "z+"],
                              "espesor_m": capas["espesor_capa_m"], "conductividad_w_m_k": capas["k"],
                              "inicio_h": max(0.0, capas["inicio_h"] - hor["activaciones_h"][i]), "inicio_desde": "colocacion",
                              "retiro_h": capas["retiros_h"][j], "retiro_desde": "inicio_obra"})
    config["clima"] = dict(db._leer_jsonc(RAIZ / db.CLIMAS[CLIMA][0])["clima"])
    config["clima"]["inicio_obra_utc"] = "2026-11-19T13:00:00Z"
    if horizonte_h > 672.0:
        config["clima"]["snapshot"] = dict(config["clima"]["snapshot"], extension="ciclo_diario")
    config["sitio"] = {"grupo": 3, "latitud": -31.1115, "longitud": -71.6147}
    o = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    plan = PlanConstruccion(PLAN, hor["activaciones_h"], hor["duraciones_h"], {i: 10.0 for i in range(1, 6)}, 1)
    vertidas = plan.aplicar(o)
    grilla = voxelizar_vertidas(vertidas, MALLAS[malla])
    return spec, hor, config, vertidas, grilla


def diagnosticos(malla, grilla, config, vertidas):
    k, rho, cp = pipe.campos_materiales(grilla, vertidas, config["materiales"])
    solver = DifusionCartesiana(forma=grilla.forma, dx_m=grilla.paso_m, conductividad_w_m_k=k, densidad_kg_m3=rho,
                                calor_especifico_j_kg_k=cp, fronteras={}, mascara_activa=(grilla.vertida_id > 0))
    fuente = ModeloFuenteTabulada.desde_source_terms(db.SOURCE_TERMS, CURVA, cola_cero=True)
    e_curva = float(np.trapezoid(fuente.qdot_w_m3, fuente.tiempos_s))
    vv = float(np.prod(grilla.paso_m))
    d = {"malla": malla, "paso_m": list(grilla.paso_m), "forma": list(grilla.forma),
         "voxeles_hormigon": int((grilla.vertida_id > 0).sum()), "volumen_voxel_m3": solver.volumen_m3,
         "areas_cara_m2": solver.area_cara_m2, "paso_estable_interior_s": solver.paso_estable_interior_s(),
         "paso_estable_conductancias_s": solver.paso_maximo_estable_s(),
         "dt_usado_max_s": min(600.0, 0.8 * solver.paso_maximo_estable_s()),
         "profundidad_S25_m": 0.025, "S25_fraccion_celda_borde": {a: 0.025 / (grilla.paso_m[i] / 2) for i, a in enumerate("xyz")},
         "E_curva_J_m3": e_curva, "por_vertida": {}}
    assert abs(solver.volumen_m3 - vv) < 1e-12
    for i in range(1, 6):
        v = float((grilla.vertida_id == i).sum() * vv)
        d["por_vertida"][i] = {"vol_IFC_m3": IFC[i], "vol_voxel_m3": v, "err_abs_m3": v - IFC[i],
                               "energia_generada_esperada_J": v * e_curva}
    d["total"] = {"vol_voxel_m3": sum(x["vol_voxel_m3"] for x in d["por_vertida"].values()),
                  "energia_generada_esperada_J": sum(x["energia_generada_esperada_J"] for x in d["por_vertida"].values())}
    return d


def correr_reducido(config, grilla, vertidas, observadores, fuente, tiempos, instantes_h):
    """Igual que demo.correr pero guardando campos completos solo en ``instantes_h`` (+ t=0 y final);
    las métricas M_T y M_delta siguen calculándose en cada paso por los observadores."""
    from mcplan.thermal.vertidas import simular_vertidas
    k, rho, cp = pipe.campos_materiales(grilla, vertidas, config["materiales"])
    final_s = float(config["tiempo_final_h"]) * 3600.0
    bordes, aire = pipe.crear_bordes_entorno(config, grilla, vertidas, final_s)
    res = simular_vertidas(
        grilla=grilla, vertidas=vertidas, conductividad_w_m_k=k, densidad_kg_m3=rho, calor_especifico_j_kg_k=cp,
        hidratacion=fuente, tiempo_final_s=final_s, fraccion_paso_estable=float(config.get("fraccion_paso_estable", 0.8)),
        paso_tiempo_max_s=float(config.get("paso_tiempo_max_s", 600.0)),
        tiempos_muestreo_s=sorted({0.0, final_s, *[h * 3600.0 for h in instantes_h]}),
        fronteras=pipe.crear_fronteras(config), seguimiento=list(observadores), bordes=bordes,
        tiempos_activacion_s_por_voxel=tiempos)
    return res, bordes, aire


def correr(malla, ai_mm=50, reducido=True, suelo_L=0.8, regla=None, horizonte_h=1344.0):
    t_ini = time.time()
    SALIDA.mkdir(parents=True, exist_ok=True)
    spec, hor, config, vertidas, grilla = preparar(malla, ai_mm, suelo_L, regla, horizonte_h)
    diag = diagnosticos(malla, grilla, config, vertidas)
    caras = ("x-", "x+", "y-", "y+", "z-", "z+")
    tiempos = dc.tiempos_frente(grilla, vertidas, hor)
    o_fin = SeguimientoMetricasDecision(grilla, 0.025, caras)
    o14 = SeguimientoMetricasDecision(grilla, 0.025, caras, tiempo_max_s=336 * 3600.0)

    class Serie:
        def __init__(self):
            self.h, self.d, self.tm, self.sig = [], [], [], 0.0

        def registrar(self, *a):
            pass

        def registrar_paso(self, t_s, dt_s, solver, T, f):
            u = o_fin.ultimo
            if u and t_s + 1e-6 >= self.sig and abs(u["t_s"] - t_s) < 1e-6:
                self.h.append(round(t_s / 3600, 3)); self.d.append(u["delta_c"]); self.tm.append(u["T_max_c"]); self.sig = t_s + 3600.0
    serie = Serie()
    fuente = ModeloFuenteTabulada.desde_source_terms(db.SOURCE_TERMS, CURVA, cola_cero=True)
    if reducido:
        res, _, _ = correr_reducido(config, grilla, vertidas, (o_fin, o14, serie), fuente, tiempos, (336.0, 860.0))
    else:
        res, _, _ = demo.correr(config, grilla, vertidas, observadores=(o_fin, o14, serie), modelo_hidratacion=fuente,
                                tiempos_activacion_s_por_voxel=tiempos)
    out = {"suelo_L": suelo_L, "regla": regla, "horizonte_h": horizonte_h,
           "serie": {"h": serie.h, "delta_c": serie.d, "Tmax_c": serie.tm}, "ai_mm": ai_mm, "muestreo": "reducido" if reducido else "horario", "diagnostico": diag, "h14": o14.resumen(), "hfin": o_fin.resumen(),
           "energia": {"fuente_J": float(res.energia_fuente_j), "activacion_J": float(res.energia_activacion_j),
                       "frontera_J": float(res.energia_frontera_j), "final_J": float(res.energia_final_j),
                       "error_balance_J": float(res.error_balance_j),
                       "fuente_vs_esperada_rel": float(res.energia_fuente_j / diag["total"]["energia_generada_esperada_J"] - 1.0)},
           "ultimo_paso_maximo_estable_s": float(res.ultimo_paso_maximo_estable_s),
           "n_muestras_observador": o_fin.n_muestras, "segundos": time.time() - t_ini}
    (SALIDA / f"malla_{malla}_AI{ai_mm}_{'red' if reducido else 'hor'}_S{suelo_L}_R{regla}_H{int(horizonte_h)}_{PLAN}_{CURVA}_{CLIMA}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("energia", "n_muestras_observador", "segundos", "ultimo_paso_maximo_estable_s")}, indent=1))
    print("MdT56", out["hfin"]["M_delta_c"], out["hfin"]["t_M_delta_h"], out["hfin"]["xyz_S25_min_m"], "Tmax", out["hfin"]["M_T_c"])
    return out


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    modo, malla = sys.argv[1:3]
    ai_mm = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    reducido = (sys.argv[4] != "horario") if len(sys.argv) > 4 else True
    suelo_L = float(sys.argv[5]) if len(sys.argv) > 5 else 0.8
    regla = (sys.argv[6] if sys.argv[6] != "none" else None) if len(sys.argv) > 6 else None
    horizonte_h = float(sys.argv[7]) if len(sys.argv) > 7 else 1344.0
    if modo == "verificar":
        _, _, config, vertidas, grilla = preparar(malla)
        print(json.dumps(diagnosticos(malla, grilla, config, vertidas), indent=1))
    else:
        correr(malla, ai_mm, reducido, suelo_L, regla, horizonte_h)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
