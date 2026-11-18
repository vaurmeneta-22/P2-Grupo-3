"""Tablas de sensibilidad (física y numérica) comparadas con la corrida de referencia.

    python scripts/defensa_sensibilidad.py
Escribe results/defensa/sensibilidad_fisica.csv y sensibilidad_numerica.csv.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "scripts"))
import defensa_figuras as fig  # noqa: E402

SALIDA = RAIZ / "results/defensa"
CURVAS_CLIMA = (("M4_SS", "fria"), ("M1_SS", "fria"))


def fila(d, ref, etiqueta, tipo):
    m, r = d["hfin"], ref["hfin"]
    sb = d.get("hfin_sin_base", {})
    return {"tipo": tipo, "variante": etiqueta, "curva": d["curva"], "clima": d["clima"],
            "D_h": d["D_h"], "dx_m": d["dx_m"], "dt_max_s": d["dt_max_s"],
            "M_T_14d_c": d["h14"]["M_T_c"], "M_delta_14d_c": d["h14"]["M_delta_c"],
            "M_T_fin_c": m["M_T_c"], "M_delta_fin_c": m["M_delta_c"], "t_M_delta_fin_h": m["t_M_delta_h"],
            "M_delta_fin_sin_base_c": sb.get("M_delta_c"),
            "dM_T_vs_ref_c": m["M_T_c"] - r["M_T_c"], "dM_delta_vs_ref_c": m["M_delta_c"] - r["M_delta_c"],
            "dM_delta_14d_vs_ref_c": d["h14"]["M_delta_c"] - ref["h14"]["M_delta_c"],
            "cumple_14d": int(d["h14"]["M_T_c"] <= 70 and d["h14"]["M_delta_c"] <= 20),
            "cumple_fin": int(m["M_T_c"] <= 70 and m["M_delta_c"] <= 20),
            "cambia_clasificacion_fin": int((m["M_delta_c"] <= 20) != (r["M_delta_c"] <= 20)),
            "M_delta_centro_voxel_fin_c": m.get("M_delta_centro_voxel_c"), "M_delta_cara_fin_c": m.get("M_delta_cara_c"),
            "balance_error_rel": abs(d["balance"]["error_J"]) / max(abs(d["balance"]["energia_fuente_J"])
                                                                    + abs(d["balance"]["energia_activacion_J"]), 1.0)}


def main() -> int:
    SALIDA.mkdir(parents=True, exist_ok=True)
    fisica, numerica = [], []
    for curva, clima in CURVAS_CLIMA:
        ref = fig.corrida("f_P2_L2m_U63", curva, clima)
        fisica.append(fila(ref, ref, "referencia f_P2_L2m_U63", "referencia"))
        for f in sorted(fig.CACHE.glob(f"s_*__{curva}__{clima}__*.json")):
            d = json.loads(f.read_text(encoding="utf-8"))
            if abs(d["horizonte_h"] - 1344.0) > 1e-6:
                continue
            fisica.append(fila(d, ref, d["plan_id"][2:], d["categoria"]))
        for pid in ("f_P2_L2m_U63", "f_P2_L0_8m_U63", "d_L2m_T10_PUR60_AI50", "d_L2m_T10_PUR80_AI0"):
            try:
                base = fig.corrida(pid, curva, clima)
            except FileNotFoundError:
                continue
            numerica.append(fila(base, base, f"{pid} (dx 0,25 m; dt 600 s)", "referencia"))
            for tag, dx, dt in (("dx0_5", 0.5, 600.0), ("dx0_125", 0.125, 600.0), ("dt300", 0.25, 300.0), ("dt150", 0.25, 150.0)):
                try:
                    d = fig.corrida(f"n_{pid}_{tag}", curva, clima, dx=dx, dt=dt)
                except FileNotFoundError:
                    continue
                numerica.append(fila(d, base, f"{pid} dx={dx:g} dt={dt:g}", "numerica"))
    for nombre, filas in (("sensibilidad_fisica", fisica), ("sensibilidad_numerica", numerica)):
        if filas:
            with (SALIDA / f"{nombre}.csv").open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(filas[0]))
                w.writeheader()
                w.writerows(filas)
            print(f"== {nombre}")
            for r in filas:
                print(f"{r['curva']:6s} {r['variante'][:46]:46s} T14={r['M_T_14d_c']:.2f} d14={r['M_delta_14d_c']:.2f} "
                      f"Tfin={r['M_T_fin_c']:.2f} dfin={r['M_delta_fin_c']:.2f}@{r['t_M_delta_fin_h']:.0f}h "
                      f"dref={r['dM_delta_vs_ref_c']:+.2f} flip={r['cambia_clasificacion_fin']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
