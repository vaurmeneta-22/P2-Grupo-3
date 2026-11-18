"""Consolida el caché de corridas en tablas (por escenario, por plan) y aplica P_f.

Interpretaciones de pesos (fijadas ANTES de comparar planes):
  nominal_M1_SS : curva M1_SS (peso 1) × 3 climas (11/31, 10/31, 10/31).
  ocho_curvas   : 8 curvas con peso uniforme 1/8 × 3 climas, independencia curva-clima.
                  Es una hipótesis de enumeración de la pauta, no se demuestra que las
                  ocho curvas sean réplicas de una misma mezcla (ver informe).
  solo_SS       : 4 curvas SS (peso 1/4) × 3 climas (campaña de mayor confianza).

    python scripts/defensa_consolidar.py --prefijo c_ --salida results/defensa/tabla_cribado
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
from mcplan.thermal.decision import ResultadoEscenario, agregar_probabilidades  # noqa: E402

CACHE = RAIZ / "results/defensa/cache"
CURVAS_SS = ("M1_SS", "M2_SS", "M3_SS", "M4_SS")
CURVAS_22A = ("M1_22A", "M2_22A", "M3_22A", "M4_22A")
INTERP = {
    "nominal_M1_SS": {"M1_SS": 1.0},
    "solo_SS": {c: 0.25 for c in CURVAS_SS},
    "ocho_curvas": {c: 0.125 for c in CURVAS_SS + CURVAS_22A},
}


def cargar(filtro) -> list[dict]:
    """Lee el caché; ante entradas repetidas (código distinto) conserva la más reciente."""
    vistas = {}
    for f in sorted(CACHE.glob("*.json"), key=lambda p: p.stat().st_mtime):
        d = json.loads(f.read_text(encoding="utf-8"))
        if filtro(d):
            vistas[(d["plan_id"], d["curva"], d["clima"], d["dx_m"], d["dt_max_s"], d["s25_m"],
                    d["horizonte_h"])] = d
    return list(vistas.values())


def filas_escenario(corridas, horizonte: str, interp: str, pesos_curva: dict) -> list[dict]:
    filas = []
    for d in corridas:
        if d["curva"] not in pesos_curva:
            continue
        m = d[horizonte]
        peso = pesos_curva[d["curva"]] * d["peso_clima"]
        filas.append({"plan_id": d["plan_id"], "escenario_id": f"{d['curva']}__{d['clima']}", "curva": d["curva"],
                      "clima": d["clima"], "peso": peso, "M_T_c": m["M_T_c"], "M_delta_c": m["M_delta_c"],
                      "t_M_T_h": m["t_M_T_h"], "t_M_delta_h": m["t_M_delta_h"],
                      "xyz_M_T_m": m["xyz_M_T_m"], "xyz_S25_min_m": m["xyz_S25_min_m"],
                      "M_delta_centro_voxel_c": m.get("M_delta_centro_voxel_c"),
                      "M_delta_cara_c": m.get("M_delta_cara_c"), "D_h": d["D_h"], "dx_m": d["dx_m"],
                      "dt_max_s": d["dt_max_s"], "horizonte": horizonte})
    return filas


def tabla_planes(corridas, interp: str, pesos_curva: dict, limites=(70.0, 20.0, 0.05)) -> list[dict]:
    salida = []
    planes = sorted({d["plan_id"] for d in corridas})
    for pid in planes:
        mias = [d for d in corridas if d["plan_id"] == pid]
        esperadas = {(c, w) for c in pesos_curva for w in ("fria", "central", "calida")}
        presentes = {(d["curva"], d["clima"]) for d in mias if d["curva"] in pesos_curva}
        if presentes != esperadas:
            continue  # escenario incompleto para esta interpretación
        fila = {"plan_id": pid, "interpretacion": interp, "categoria": mias[0].get("categoria", ""),
                "D_h": mias[0]["D_h"], "n_escenarios": len(esperadas)}
        hs = [("h14", "14d"), ("h28", "28d")]
        if all("hfin" in d and d["horizonte_h"] > 672.0 for d in mias):
            hs.append(("hfin", "fin"))
        for horizonte, etiqueta in hs:
            esc = filas_escenario(mias, horizonte, interp, pesos_curva)
            res = [ResultadoEscenario(f["plan_id"], f["escenario_id"], f["peso"], f["M_T_c"], f["M_delta_c"],
                                      f["t_M_T_h"], f["t_M_delta_h"]) for f in esc]
            agg = agregar_probabilidades(res, limite_T_c=limites[0], limite_delta_c=limites[1],
                                         limite_P_f=limites[2])[0]
            fila.update({f"P_T_{etiqueta}": agg["P_T"], f"P_delta_{etiqueta}": agg["P_delta"],
                         f"P_f_{etiqueta}": agg["P_f"], f"cumple_{etiqueta}": agg["cumple_P_f"],
                         f"max_M_T_{etiqueta}": max(f["M_T_c"] for f in esc),
                         f"max_M_delta_{etiqueta}": max(f["M_delta_c"] for f in esc),
                         f"n_fallos_{etiqueta}": sum(1 for f in esc if f["M_T_c"] > limites[0] or f["M_delta_c"] > limites[1])})
        fila["cumple_ambos_horizontes"] = bool(fila["cumple_14d"] and fila["cumple_28d"]
                                               and fila.get("cumple_fin", True))
        salida.append(fila)
    return salida


def escenarios_anchos(corr) -> list[dict]:
    """Una fila por plan × curva × clima con ambos horizontes, instantes, ubicaciones y pesos."""
    filas = []
    for d in sorted(corr, key=lambda x: (x["plan_id"], x["curva"], x["clima"])):
        fila = {"plan_id": d["plan_id"], "curva": d["curva"], "clima": d["clima"], "D_h": d["D_h"],
                "peso_clima": d["peso_clima"],
                "peso_nominal_M1_SS": d["peso_clima"] if d["curva"] == "M1_SS" else 0.0,
                "peso_ocho_curvas": 0.125 * d["peso_clima"],
                "peso_solo_SS": 0.25 * d["peso_clima"] if d["curva"] in CURVAS_SS else 0.0,
                "horizonte_h": d["horizonte_h"], "dx_m": d["dx_m"], "dt_max_s": d["dt_max_s"]}
        for hz, et in (("h14", "14d"), ("h28", "28d"), ("hfin", "fin")):
            if hz not in d or (hz == "hfin" and d["horizonte_h"] <= 672.0):
                continue
            m = d[hz]
            xyz_t, xyz_s = m["xyz_M_T_m"], m["xyz_S25_min_m"]
            fila.update({
                f"M_T_{et}_c": m["M_T_c"], f"t_M_T_{et}_h": m["t_M_T_h"], f"x_M_T_{et}_m": xyz_t[0],
                f"y_M_T_{et}_m": xyz_t[1], f"z_M_T_{et}_m": xyz_t[2],
                f"M_delta_{et}_c": m["M_delta_c"], f"t_M_delta_{et}_h": m["t_M_delta_h"],
                f"x_S25min_{et}_m": xyz_s[0], f"y_S25min_{et}_m": xyz_s[1], f"z_S25min_{et}_m": xyz_s[2],
                f"I_T_{et}": int(m["M_T_c"] > 70.0), f"I_delta_{et}": int(m["M_delta_c"] > 20.0),
                f"I_f_{et}": int(m["M_T_c"] > 70.0 or m["M_delta_c"] > 20.0),
                f"M_delta_centro_voxel_{et}_c": m.get("M_delta_centro_voxel_c"),
                f"M_delta_cara_{et}_c": m.get("M_delta_cara_c")})
        b = d["balance"]
        fila["balance_error_rel"] = abs(b["error_J"]) / max(abs(b["energia_fuente_J"]) + abs(b["energia_activacion_J"]), 1.0)
        fila["balance_error_J"] = b["error_J"]
        filas.append(fila)
    return filas


def escribir(ruta: Path, filas: list[dict]):
    if not filas:
        return
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefijo", default="")
    ap.add_argument("--ids", default="")
    ap.add_argument("--salida", required=True)
    ap.add_argument("--dx", type=float, default=0.25)
    ap.add_argument("--dt", type=float, default=600.0)
    a = ap.parse_args()
    ids = set(a.ids.split(",")) if a.ids else None

    def filtro(d):
        return (d["plan_id"].startswith(a.prefijo) and (ids is None or d["plan_id"] in ids)
                and abs(d["dx_m"] - a.dx) < 1e-9 and abs(d["dt_max_s"] - a.dt) < 1e-9)

    corr = cargar(filtro)
    salida = Path(a.salida) if Path(a.salida).is_absolute() else RAIZ / a.salida
    todas = []
    for interp, pesos in INTERP.items():
        todas += tabla_planes(corr, interp, pesos)
    escribir(salida / "comparacion_planes.csv", todas)
    esc = escenarios_anchos(corr)
    escribir(salida / "resultados_por_escenario.csv", esc)
    print(f"{len(corr)} corridas -> {salida}")
    for f in todas:
        if f["interpretacion"] == "nominal_M1_SS":
            print(f"{f['plan_id']}: D={f['D_h']:.1f} h  P_f14={f['P_f_14d']:.3f} P_f28={f['P_f_28d']:.3f}  "
                  f"maxΔ14={f['max_M_delta_14d']:.2f} maxΔ28={f['max_M_delta_28d']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
