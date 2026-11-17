"""Prepara tres trayectorias GEFS completas y trazables para la evaluación final.

Selecciona un miembro real representativo de cada tercil de temperatura media
durante la ventana de 14 días. Cada miembro conserva toda su secuencia temporal;
los pesos son las fracciones de miembros de cada tercil, no pesos arbitrarios.
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import demo_pipeline_ifc_3d as pipe
from mcplan.weather.serie import a_utc, iso_utc

ENSEMBLE = ROOT / "data/weather/2026-11-16T213220Z_GEFS_ensemble"
DOSSIER = ROOT / "configs/evaluacion_final.jsonc"
SALIDA_CLIMA = ROOT / "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias"
CONFIG_CLIMAS = ROOT / "configs/climas/gefs_20261116"
T0_LOCAL = "2026-11-19T10:00:00-03:00"
T0_UTC = "2026-11-19T13:00:00Z"
CLIMA_IDS = {
    "fria": "GEFS_20261116_frio",
    "central": "GEFS_20261116_central",
    "calida": "GEFS_20261116_calido",
}


def leer_jsonc(path: Path) -> dict:
    return json.loads(pipe.quitar_comentarios_jsonc(path.read_text(encoding="utf-8")))


def seleccionar_grupos(miembros: np.ndarray) -> list[dict]:
    medias = np.mean(miembros, axis=1)
    orden = np.argsort(medias, kind="stable")
    etiquetas = ("fria", "central", "calida")
    grupos = np.array_split(orden, 3)
    mediana_ensamble = np.median(miembros, axis=0)
    seleccion = []
    for etiqueta, indices in zip(etiquetas, grupos):
        centro = np.median(miembros[indices], axis=0)
        if etiqueta == "central":
            centro = mediana_ensamble
        distancias = np.sqrt(np.mean((miembros[indices] - centro) ** 2, axis=1))
        elegido = int(indices[int(np.argmin(distancias))])
        seleccion.append({
            "etiqueta": etiqueta,
            "indices": [int(i) for i in indices],
            "indice_elegido": elegido,
            "peso": float(len(indices) / len(miembros)),
            "media_grupo_min_c": float(np.min(medias[indices])),
            "media_grupo_max_c": float(np.max(medias[indices])),
            "media_representante_c": float(medias[elegido]),
            "rmse_hacia_centro_c": float(np.min(distancias)),
        })
    return seleccion


def guardar_snapshot(item: dict, keys: list[str], horas: list[datetime],
                     valores: np.ndarray, digest_original: str,
                     ensemble_metadata: dict) -> tuple[Path, dict]:
    etiqueta = item["etiqueta"]
    carpeta = SALIDA_CLIMA / etiqueta
    carpeta.mkdir(parents=True, exist_ok=True)
    member_id = keys[item["indice_elegido"]]
    member_label = "control" if member_id == "temperature_2m" else member_id.removeprefix("temperature_2m_")
    forecast = {
        "time_utc": [iso_utc(t) for t in horas],
        "temperature_2m": [round(float(v), 4) for v in valores[item["indice_elegido"]]],
        "units": {"temperature_2m": "degC"},
    }
    (carpeta / "forecast.json").write_text(
        json.dumps(forecast, indent=2), encoding="utf-8")

    grupo_ids = ["control" if keys[i] == "temperature_2m"
                 else keys[i].removeprefix("temperature_2m_") for i in item["indices"]]
    metadata = {
        "source": "Open-Meteo Ensemble API; NOAA Global Ensemble Forecast System",
        "model": ensemble_metadata["model"],
        "member": member_label,
        "selection": f"miembro medoid del tercil {etiqueta} por temperatura media de la ventana",
        "represented_member_ids": grupo_ids,
        "represented_member_count": len(grupo_ids),
        "weight": item["peso"],
        "group_mean_temperature_range_c": [item["media_grupo_min_c"], item["media_grupo_max_c"]],
        "selected_member_mean_temperature_c": item["media_representante_c"],
        "rmse_to_group_or_ensemble_centre_c": item["rmse_hacia_centro_c"],
        "ensemble_member_count": len(keys),
        "site": "Grupo 3: Mantos de Hornillos / Ovalle, Camarico – Ecos del Mar, Coquimbo",
        "latitude": ensemble_metadata["coordinates"]["latitude"],
        "longitude": ensemble_metadata["coordinates"]["longitude"],
        "timezone": "UTC",
        "coordinates_source": "Tabla docente de sitios del curso; mismos datos usados en E11",
        "issued_at": None,
        "issued_at_note": "la API no expone la hora de emisión del ciclo GEFS",
        "retrieved_at": ensemble_metadata["retrieved_at"],
        "valid_from": iso_utc(horas[0]),
        "valid_to": iso_utc(horas[-1]),
        "variables": ["temperature_2m"],
        "units": {"temperature_2m": "degC"},
        "original_file": "../../2026-11-16T213220Z_GEFS_ensemble/original.json",
        "original_sha256": digest_original,
        "ensemble_manifest_file": "../../2026-11-16T213220Z_GEFS_ensemble/ensemble_metadata.json",
        "series_kind": "pronostico",
        "synthetic": False,
        "native_temporal_resolution": "3-hourly GEFS output; hourly time series interpolated by Open-Meteo",
    }
    (carpeta / "metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    climate_id = CLIMA_IDS[etiqueta]
    config_rel = f"configs/climas/gefs_20261116/{etiqueta}.jsonc"
    climate_config = {
        "clima": {
            "fuente": "snapshot",
            "inicio_obra_utc": T0_UTC,
            "snapshot": {
                "carpeta": carpeta.relative_to(ROOT).as_posix(),
                "raiz": "data/weather",
                "extension": None,
            },
        }
    }
    CONFIG_CLIMAS.mkdir(parents=True, exist_ok=True)
    (ROOT / config_rel).write_text(
        json.dumps(climate_config, indent=2) + "\n", encoding="utf-8")
    return carpeta, {"id": climate_id, "config_jsonc": config_rel, "peso": item["peso"]}


def actualizar_dossier(climas: list[dict], ids: list[str]) -> None:
    texto = DOSSIER.read_text(encoding="utf-8")
    texto, n_t0 = re.subn(
        r'"inicio_obra_utc"\s*:\s*"[^"]+"',
        f'"inicio_obra_utc": "{T0_UTC}"', texto, count=1)
    if n_t0 != 1:
        raise ValueError("No se encontró inicio_obra_utc para actualizar")
    texto, n_local = re.subn(
        r'("inicio_obra_local"\s*:\s*)"[^"]+"',
        rf'\g<1>"{T0_LOCAL}"', texto, count=1)
    if n_local == 0:
        texto = texto.replace(
            f'"inicio_obra_utc": "{T0_UTC}",',
            f'"inicio_obra_utc": "{T0_UTC}",\n  "inicio_obra_local": "{T0_LOCAL}", // 10:00 hora de Ovalle (UTC-3).')
    # Fija coordenadas iguales a las usadas en el archivo meteorológico.
    texto = re.sub(r'"latitud"\s*:\s*null', '"latitud": -31.1115', texto, count=1)
    texto = re.sub(r'"longitud"\s*:\s*null', '"longitud": -71.6147', texto, count=1)
    texto = texto.replace('"fuente_coordenadas": "Coordenadas pendientes de sitios.yaml"',
                          '"fuente_coordenadas": "Tabla docente de sitios del curso; ver metadata GEFS"')
    block = '"climas": [\n'
    for c, etiqueta in zip(climas, ("fría", "central", "cálida")):
        block += (f'    // Trayectoria GEFS {etiqueta}; miembro completo representativo de un tercil.\n'
                  f'    {{"id":"{c["id"]}", "config_jsonc":"{c["config_jsonc"]}", "peso":{c["peso"]:.12f}}},\n')
    block = block.rstrip(",\n") + '\n  ]'
    texto, n_climas = re.subn(
        r'"climas"\s*:\s*\[[\s\S]*?\n  \](?=,\n  "pesos_curva_clima_independientes")',
        block, texto, count=1)
    if n_climas != 1:
        raise ValueError("No se localizó el bloque climas del dossier")
    central_id = ids[1]
    texto = re.sub(r'"escenario_id"\s*:\s*"[^"]+"',
                   f'"escenario_id": "M1_SS__{central_id}"', texto, count=1)
    DOSSIER.write_text(texto, encoding="utf-8")


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    datos = json.loads((ENSEMBLE / "original.json").read_text(encoding="utf-8"))
    info = json.loads((ENSEMBLE / "ensemble_metadata.json").read_text(encoding="utf-8"))
    hourly = datos["hourly"]
    all_times = [a_utc(t) for t in hourly["time"]]
    t0 = a_utc(T0_UTC)
    tf = t0 + timedelta(hours=336.0)
    indices = [i for i, t in enumerate(all_times) if t0 <= t <= tf]
    expected = 337
    if len(indices) != expected or all_times[indices[0]] != t0 or all_times[indices[-1]] != tf:
        raise ValueError("El ensamble no cubre exactamente las 337 horas de t0 a t0+14 días")

    keys = ["temperature_2m"] + [f"temperature_2m_member{i:02d}" for i in range(1, 31)]
    matrix_all = np.asarray([
        [np.nan if v is None else float(v) for v in hourly[key]] for key in keys
    ], dtype=float)
    matrix = matrix_all[:, indices]
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Hay miembros/horas faltantes en la ventana climática evaluada")
    grupos = seleccionar_grupos(matrix)
    digest = hashlib.sha256((ENSEMBLE / "original.json").read_bytes()).hexdigest()

    SALIDA_CLIMA.mkdir(parents=True, exist_ok=True)
    snapshots, climate_configs = [], []
    for grupo in grupos:
        carpeta, climate_cfg = guardar_snapshot(
            grupo, keys, [all_times[i] for i in indices], matrix, digest, info)
        snapshots.append(carpeta)
        climate_configs.append(climate_cfg)

    # Manifiesto de reducción: miembros representados y extremos que quedan fuera.
    medias = np.mean(matrix, axis=1)
    ranking = np.argsort(medias)
    informe = {
        "method": "rank 31 complete trajectories by mean temperature over the same 14-day window; split into three terciles; select each tercile medoid, central medoid closest to pointwise ensemble median",
        "t0_local": T0_LOCAL,
        "t0_utc": T0_UTC,
        "window_end_utc": iso_utc(tf),
        "member_count": len(keys),
        "weights": {g["etiqueta"]: g["peso"] for g in grupos},
        "selected": [
            {"group": g["etiqueta"], "member": keys[g["indice_elegido"]],
             "weight": g["peso"],
             "represented_member_ids": [keys[i] for i in g["indices"]],
             "group_mean_range_c": [g["media_grupo_min_c"], g["media_grupo_max_c"]],
             "selected_mean_c": g["media_representante_c"]}
            for g in grupos
        ],
        "full_ensemble_mean_extremes_c": {
            "coldest_member": keys[int(ranking[0])], "coldest_mean_c": float(medias[ranking[0]]),
            "warmest_member": keys[int(ranking[-1])], "warmest_mean_c": float(medias[ranking[-1]]),
        },
        "notice": "The three selected medoids do not retain every member extreme; the excluded full-ensemble min/max are reported for sensitivity.",
    }
    (SALIDA_CLIMA / "selection_metadata.json").write_text(
        json.dumps(informe, indent=2, ensure_ascii=False), encoding="utf-8")
    actualizar_dossier(climate_configs, [c["id"] for c in climate_configs])
    print(f"Inicio local/UTC: {T0_LOCAL} / {T0_UTC}")
    print(f"Ventana: {T0_UTC} a {iso_utc(tf)} ({expected} registros horarios)")
    for g, cfg in zip(grupos, climate_configs):
        print(f"{g['etiqueta']}: member={keys[g['indice_elegido']]}, "
              f"peso={g['peso']:.6f}, promedio={g['media_representante_c']:.2f} °C")
    print(f"Configuración actualizada: {DOSSIER}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
