"""Audita las ocho series Qdot de data/hidratacion/source_terms.txt.

Comprueba unidades, origen temporal, cobertura, calor total frente al encabezado del
archivo (beta = 2,5 W/m²K) y cola truncada. Escribe results/defensa/curvas_auditoria.csv.
"""
from __future__ import annotations

import csv
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcplan.thermal.hidratacion import ModeloFuenteTabulada  # noqa: E402

ARCHIVO = ROOT / "data/hidratacion/source_terms.txt"
SALIDA = ROOT / "results/defensa"
RHO_CP_ENSAYO = 2476.0 * 1023.0   # del encabezado del archivo (probeta)
RHO_CP_MODELO = 2400.0 * 900.0    # configs/semana3_fundacion.jsonc (M1)


def encabezado_totales() -> dict:
    out = {}
    for linea in ARCHIVO.read_text(encoding="utf-8").splitlines():
        m = re.match(r"#\s+(M\d_(?:SS|22A))\s+([\d.]+) / ([\d.]+) / ([\d.]+)\s+([\d.]+) / ([\d.]+) / ([\d.]+)", linea)
        if m:
            out[m.group(1)] = {"MJ_m3_beta2_5": float(m.group(5)), "J_g_beta2_5": float(m.group(6)),
                               "K_adiabatico_beta2_5": float(m.group(7))}
    return out


def main() -> int:
    SALIDA.mkdir(parents=True, exist_ok=True)
    cab = encabezado_totales()
    filas = []
    for col in ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"):
        m = ModeloFuenteTabulada.desde_source_terms(ARCHIVO, col, cola_cero=True)
        t_h = np.asarray(m.tiempos_s) / 3600.0
        q = np.asarray(m.qdot_w_m3)
        e = float(np.trapezoid(q, np.asarray(m.tiempos_s))) / 1e6   # MJ/m3
        ultimas24 = t_h >= t_h[-1] - 24.0
        e24 = float(np.trapezoid(q[ultimas24], np.asarray(m.tiempos_s)[ultimas24])) / 1e6
        ref = cab.get(col, {})
        filas.append({
            "curva": col, "t_ini_h": t_h[0], "t_fin_h": t_h[-1], "n_muestras": len(t_h),
            "paso_h": float(np.median(np.diff(t_h))), "Q_inicial_W_m3": q[0], "Q_pico_W_m3": q.max(),
            "t_pico_h": float(t_h[q.argmax()]), "Q_final_W_m3": q[-1],
            "calor_total_MJ_m3": e, "calor_encabezado_MJ_m3": ref.get("MJ_m3_beta2_5"),
            "dif_vs_encabezado_pct": (100 * (e / ref["MJ_m3_beta2_5"] - 1)) if ref else None,
            "calor_g_cemento_J_g": e * 1e6 / 350.77 / 1000.0,
            "ascenso_adiabatico_ensayo_K": e * 1e6 / RHO_CP_ENSAYO,
            "ascenso_adiabatico_modelo_K": e * 1e6 / RHO_CP_MODELO,
            "calor_ultimas_24h_MJ_m3": e24, "fraccion_ultimas_24h_pct": 100 * e24 / e,
            "q_final_sobre_pico_pct": 100 * q[-1] / q.max(),
        })
    with (SALIDA / "curvas_auditoria.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    for r in filas:
        print(f"{r['curva']}: fin={r['t_fin_h']:.1f} h, pico={r['Q_pico_W_m3']:.0f} W/m3 @ {r['t_pico_h']:.1f} h, "
              f"E={r['calor_total_MJ_m3']:.2f} MJ/m3 (encab {r['calor_encabezado_MJ_m3']}; dif {r['dif_vs_encabezado_pct']:.2f}%), "
              f"dTad modelo={r['ascenso_adiabatico_modelo_K']:.1f} K, Q_fin/pico={r['q_final_sobre_pico_pct']:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
