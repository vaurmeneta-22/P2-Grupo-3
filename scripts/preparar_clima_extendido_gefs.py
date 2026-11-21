"""Prepara las tres trayectorias GEFS seleccionadas a 28 días o al fin común del archivo.

Usa ``original.json`` (descarga del 16-11-2026, 35 días de pronóstico) y los MISMOS
miembros representativos (frío/central/cálido) y pesos de
``selection_metadata.json``; no se re-seleccionan miembros ni se modifican pesos.
Las primeras 337 h coinciden con las trayectorias de 14 días (se verifica).
Más allá de ~15 días el pronóstico GEFS tiene poca habilidad: esa cola es una
trayectoria del ensamble, no un pronóstico fiable, y se etiqueta así.

    python scripts/preparar_clima_extendido_gefs.py
    python scripts/preparar_clima_extendido_gefs.py --archivo-completo
"""
from __future__ import annotations

from datetime import timedelta
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcplan.weather.serie import a_utc, iso_utc  # noqa: E402

ENS = ROOT / "data/weather/2026-11-16T213220Z_GEFS_ensemble"
BASE = ROOT / "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias"
SALIDA = ROOT / "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias_28d"
CFG = ROOT / "configs/climas/gefs_20261116_28d"
T0 = "2026-11-19T13:00:00Z"
HORAS = 672


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archivo-completo", action="store_true")
    args = parser.parse_args()
    salida, cfg, horas = SALIDA, CFG, HORAS
    datos = json.loads((ENS / "original.json").read_text(encoding="utf-8"))
    hourly = datos["hourly"]
    tiempos = [a_utc(t) for t in hourly["time"]]
    t0 = a_utc(T0)
    sel = json.loads((BASE / "selection_metadata.json").read_text(encoding="utf-8"))
    if args.archivo_completo:
        miembros = [item["member"] for item in sel["selected"]]
        fin_comun = next((t for i, t in enumerate(tiempos) if t >= t0
                          and any(hourly[m][i] is None for m in miembros)), None)
        fin_comun = tiempos[-1] if fin_comun is None else fin_comun - timedelta(hours=1)
        horas = int((fin_comun - t0).total_seconds() / 3600)
        salida = BASE.with_name(BASE.name + "_archivo")
        cfg = CFG.with_name("gefs_20261116_archivo")
    tf = t0 + timedelta(hours=horas)
    idx = [i for i, t in enumerate(tiempos) if t0 <= t <= tf]
    if len(idx) != horas + 1 or tiempos[idx[-1]] != tf:
        raise ValueError("original.json no cubre el intervalo solicitado")
    digest = hashlib.sha256((ENS / "original.json").read_bytes()).hexdigest()
    salida.mkdir(parents=True, exist_ok=True)
    cfg.mkdir(parents=True, exist_ok=True)
    resumen = {}
    for item in sel["selected"]:
        etiqueta, clave = item["group"], item["member"]
        vals = [None if hourly[clave][i] is None else float(hourly[clave][i]) for i in idx]
        if any(v is None for v in vals):
            raise ValueError(f"{clave}: faltan horas en la ventana solicitada")
        previo = json.loads((BASE / etiqueta / "forecast.json").read_text(encoding="utf-8"))
        dif = np.max(np.abs(np.asarray(previo["temperature_2m"]) - np.asarray(vals[:337])))
        if dif > 1e-3:
            raise ValueError(f"{etiqueta}: las primeras 337 h difieren de la trayectoria de 14 días ({dif})")
        if args.archivo_completo:
            anterior = json.loads((SALIDA / etiqueta / "forecast.json").read_text())
            if anterior["temperature_2m"] != [round(v, 4) for v in vals[:673]]:
                raise ValueError(f"{etiqueta}: cambió la ventana original de 672 h")
        carpeta = salida / etiqueta
        carpeta.mkdir(exist_ok=True)
        (carpeta / "forecast.json").write_text(json.dumps({
            "time_utc": [iso_utc(tiempos[i]) for i in idx],
            "temperature_2m": [round(v, 4) for v in vals],
            "units": {"temperature_2m": "degC"}}, indent=2), encoding="utf-8")
        meta = json.loads((BASE / etiqueta / "metadata.json").read_text(encoding="utf-8"))
        meta.update({
            "valid_to": iso_utc(tf),
            "selection": meta["selection"] + "; miembro y pesos fijados con la ventana de 14 días y mantenidos sin cambio en la ventana de 28 días",
            "horizon_note": "Horas 0-336: pronóstico GEFS de los 14 días evaluados. Horas 336-672: misma trayectoria del miembro GEFS más allá de la habilidad útil de pronóstico (>~15 d); se usa solo para comprobar el comportamiento tardío, no como predicción.",
            "original_sha256": digest,
            "derived_from": "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias (primeras 337 h verificadas idénticas)",
        })
        if args.archivo_completo:
            meta["horizon_note"] = (
                f"Horas 0-336: ventana original. Horas 336-{horas}: cola archivada "
                "del mismo miembro GEFS, usada como escenario de estrés."
            )
            meta["selection"] = BASE.name + ": mismos miembros seleccionados y pesos"
            meta["derived_from"] = ENS.name + "/original.json; primeras 673 muestras verificadas"
        (carpeta / "metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
        (cfg / f"{etiqueta}.jsonc").write_text(json.dumps({"clima": {
            "fuente": "snapshot", "inicio_obra_utc": T0,
            "snapshot": {"carpeta": carpeta.relative_to(ROOT).as_posix(), "raiz": "data/weather",
                         "extension": None}}}, indent=2) + "\n", encoding="utf-8")
        resumen[etiqueta] = {"miembro": clave, "media_14d_c": float(np.mean(vals[:337])),
                             "media_15_28d_c": float(np.mean(vals[337:673])),
                             "min_c": float(np.min(vals)), "max_c": float(np.max(vals))}
    nombre_resumen = "selection_metadata_archivo.json" if args.archivo_completo else "selection_metadata_28d.json"
    (salida / nombre_resumen).write_text(json.dumps(
        {"t0_utc": T0, "horas": horas, "miembros_y_pesos": "idénticos a selection_metadata.json",
         "resumen": resumen}, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(resumen, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
