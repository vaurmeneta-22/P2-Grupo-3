"""Evalúa planes constructivos con curvas de hidratación y climas ponderados.

Lee ``configs/evaluacion_final.jsonc``. Los datos de curva y clima se congelan en
archivos locales; no consulta servicios en línea. Cada corrida reutiliza el solver
3D IFC existente y registra M_T/M_Δ en cada paso, no solo en los estados horarios.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(RAIZ / "src"))

import demo_e10_entorno as demo
import demo_pipeline_ifc_3d as pipe
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import (
    PlanConstruccion, ResultadoEscenario, SeguimientoMetricasDecision,
    agregar_probabilidades, guardar_resultados_escenario_csv,
)
from mcplan.thermal.hidratacion import ModeloFuenteTabulada
from mcplan.weather.loader import construir_aire


def leer_jsonc(ruta: Path) -> dict:
    return json.loads(pipe.quitar_comentarios_jsonc(ruta.read_text(encoding="utf-8")))


def cargar_configuracion_dossier(ruta: Path) -> tuple[dict, dict]:
    dossier = leer_jsonc(ruta)
    capas = dossier.get("config_base", ["configs/semana3_fundacion.jsonc", "configs/e10_demo.jsonc"])
    config = {}
    for capa in capas:
        config.update(leer_jsonc((RAIZ / capa).resolve()))
    # Los parámetros del dossier (horizonte, malla, paso de tiempo) prevalecen.
    config.update({k: v for k, v in dossier.items() if k not in ("config_base", "planes", "curvas", "climas", "dossier")})
    return config, dossier


def cargar_plan(cfg: dict, vertidas) -> PlanConstruccion:
    ids = {int(v.vertida_id) for v in vertidas}
    if cfg.get("usar_horario_ifc", False):
        activaciones = {int(v.vertida_id): float(v.inicio_h) for v in vertidas}
        duraciones = {int(v.vertida_id): 0.0 for v in vertidas}
        nota = cfg.get("nota_duracion", "")
        if not nota:
            raise ValueError(f"Plan {cfg.get('id')}: documenta la aproximación de duración instantánea")
    else:
        activaciones = {int(i): float(t) for i, t in (cfg.get("activaciones_h") or {}).items()}
        duraciones = {int(i): float(t) for i, t in (cfg.get("duraciones_h") or {}).items()}
    if set(activaciones) != ids or set(duraciones) != ids:
        raise ValueError(f"Plan {cfg.get('id')}: define activación y duración para IDs {sorted(ids)}")
    tcol = {int(i): float(t) for i, t in (cfg.get("temperaturas_colocacion_c") or {}).items()}
    return PlanConstruccion(str(cfg["id"]), activaciones, duraciones, tcol,
                            cfg.get("frentes_simultaneos_max"))


def leer_curvas(dossier: dict, campo: str = "curvas") -> list[tuple[dict, ModeloFuenteTabulada]]:
    curvas = []
    for cfg_original in dossier.get(campo, []):
        cfg = dict(cfg_original)
        # La corrida de sensibilidad es condicional y se guarda aparte; no
        # participa en las probabilidades nominales.
        cfg.setdefault("peso", 1.0)
        ruta = (RAIZ / cfg["archivo_source_terms"]).resolve()
        curvas.append((cfg, ModeloFuenteTabulada.desde_source_terms(
            ruta, cfg["columna"], cola_cero=bool(cfg.get("cola_cero", True)))))
    return curvas


def comprobar_pesos(items, nombre):
    if not items:
        raise ValueError(f"Faltan {nombre}; agrega sus archivos y pesos en configs/evaluacion_final.jsonc")
    ids = [str(x.get("id", "")) for x in items]
    if any(not i for i in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"Los identificadores de {nombre} deben ser únicos y no vacíos")
    pesos = [float(x["peso"]) for x in items]
    if any(not np.isfinite(p) or p <= 0 for p in pesos) or not np.isclose(sum(pesos), 1.0, atol=1e-9, rtol=0):
        raise ValueError(f"Los pesos de {nombre} deben ser positivos y sumar 1")


def construir_matriz(curvas, climas, dossier):
    """Arma pares conjuntos explícitos o producto independiente declarado."""
    por_curva = {str(c["id"]): (c, modelo) for c, modelo in curvas}
    por_clima = {str(c["id"]): c for c in climas}
    conjuntos = dossier.get("escenarios_conjuntos", [])
    if conjuntos:
        ids = [str(x.get("id", "")) for x in conjuntos]
        pesos = [float(x["peso"]) for x in conjuntos]
        if not all(ids) or len(ids) != len(set(ids)):
            raise ValueError("Los escenarios conjuntos requieren IDs únicos")
        if any(not np.isfinite(p) or p <= 0 for p in pesos) or not np.isclose(sum(pesos), 1.0, atol=1e-9, rtol=0):
            raise ValueError("Los pesos de escenarios conjuntos deben ser positivos y sumar 1")
        matriz = []
        for item in conjuntos:
            curva_id, clima_id = str(item["curva_id"]), str(item["clima_id"])
            if curva_id not in por_curva or clima_id not in por_clima:
                raise ValueError(f"Escenario {item['id']}: curva o clima inexistente")
            matriz.append((str(item["id"]), por_curva[curva_id], por_clima[clima_id], float(item["peso"])))
        return matriz
    if not dossier.get("pesos_curva_clima_independientes", False):
        raise ValueError("Define escenarios_conjuntos con pesos conjuntos, o declara explícitamente "
                         "pesos_curva_clima_independientes=true")
    comprobar_pesos([c for c, _ in curvas], "curvas de fuente volumétrica")
    comprobar_pesos(climas, "trayectorias climáticas")
    return [(f"{cc['id']}__{wx['id']}", (cc, modelo), wx,
             float(cc["peso"]) * float(wx["peso"]))
            for cc, modelo in curvas for wx in climas]


def sha256(ruta: Path) -> str:
    h = hashlib.sha256()
    with ruta.open("rb") as f:
        for bloque in iter(lambda: f.read(1024 * 1024), b""):
            h.update(bloque)
    return h.hexdigest()


def main() -> int:
    # PowerShell suele usar una página de códigos que no admite algunos símbolos
    # del resumen final (por ejemplo, ≤). Evita que una corrida completa falle
    # después de haber guardado sus resultados.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    sensibilidad = False
    salida_override = None
    argumento_config = None
    argumentos = iter(sys.argv[1:])
    for argumento in argumentos:
        if argumento == "--sensibilidad":
            sensibilidad = True
        elif argumento == "--salida":
            try:
                salida_override = next(argumentos)
            except StopIteration as exc:
                raise ValueError("--salida requiere una ruta de carpeta") from exc
        elif argumento_config is None:
            argumento_config = argumento
        else:
            raise ValueError("Uso: evaluar_dossier_final.py [config.jsonc] [--sensibilidad] [--salida carpeta]")
    ruta_cfg = Path(argumento_config) if argumento_config else RAIZ / "configs/evaluacion_final.jsonc"
    if not ruta_cfg.is_absolute():
        ruta_cfg = RAIZ / ruta_cfg
    config, dossier = cargar_configuracion_dossier(ruta_cfg)
    curvas = leer_curvas(dossier, "curvas_sensibilidad" if sensibilidad else "curvas")
    climas = dossier.get("climas", [])
    planes_cfg = dossier.get("planes", [])
    if not curvas:
        raise ValueError("No hay curvas de fuente volumétrica configuradas. Revisa el archivo source_terms del profesor.")
    if not sensibilidad and dossier.get("exigir_ocho_curvas", True) and len(curvas) != 8:
        raise ValueError(f"La evaluación final requiere las 8 series originales de Qdot; se encontraron {len(curvas)}")
    if not climas:
        raise ValueError("No hay trayectorias climáticas; agrega series completas para los 14 días.")
    matriz = construir_matriz(curvas, climas, dossier)
    if not planes_cfg:
        raise ValueError("No hay planes configurados; agrega el propuesto y las alternativas constructibles.")

    limite = dossier.get("limites", {})
    limite_t = float(limite.get("temperatura_maxima_c", 70.0))
    limite_delta = float(limite.get("diferencia_interior_superficie_c", 20.0))
    limite_pf = float(limite.get("probabilidad_fallo_maxima", 0.05))
    t0 = dossier.get("inicio_obra_utc")
    if not t0 or not dossier.get("zona_horaria_local"):
        raise ValueError("Declara inicio_obra_utc (t0 común) y zona_horaria_local")
    config["tiempo_final_h"] = float(dossier.get("tiempo_final_h", 336.0))
    config["paso_voxel_m"] = float(dossier.get("paso_voxel_m", 0.25))
    config["paso_tiempo_max_s"] = float(dossier.get("paso_tiempo_max_s", 600.0))
    climas_config = {}
    for clima_cfg in climas:
        overlay = leer_jsonc((RAIZ / clima_cfg["config_jsonc"]).resolve())
        clima_local = dict(overlay.get("clima") or {})
        if clima_local.get("fuente") != "snapshot":
            raise ValueError(f"Clima {clima_cfg['id']}: config debe usar fuente 'snapshot'")
        if (clima_local.get("snapshot") or {}).get("extension") not in (None, ""):
            raise ValueError(f"Clima {clima_cfg['id']}: no se permite extender cobertura en la corrida final")
        clima_local["inicio_obra_utc"] = t0
        aire_prueba = construir_aire(clima_local,
                                     tiempo_final_s=config["tiempo_final_h"] * 3600.0,
                                     raiz_repo=RAIZ)
        if aire_prueba.es_sintetica:
            raise ValueError(f"Clima {clima_cfg['id']}: el snapshot está marcado como sintético")
        climas_config[str(clima_cfg["id"])] = overlay

    ruta_ifc = (RAIZ / config["ifc"]).resolve()
    originales = extraer_vertidas(ruta_ifc)
    planes = [cargar_plan(p, originales) for p in planes_cfg]
    if len({p.plan_id for p in planes}) != len(planes):
        raise ValueError("Los identificadores de planes deben ser únicos")
    categorias = {str(p["id"]): str(p.get("categoria", "")) for p in planes_cfg}
    requeridas = set(dossier.get("categorias_plan_requeridas", []))
    presentes = set(categorias.values())
    if not requeridas <= presentes:
        faltan = ", ".join(sorted(requeridas - presentes))
        raise ValueError(f"Faltan categorías de planes requeridas para la comparación: {faltan}")
    resultados = []
    salida_cfg = Path(salida_override) if salida_override else Path(dossier.get("salida", "results/dossier_final"))
    salida_base = (salida_cfg if salida_cfg.is_absolute() else RAIZ / salida_cfg).resolve()
    salida = salida_base / "sensibilidad_M1_22A" if sensibilidad else salida_base
    if sensibilidad and len(curvas) != 1:
        raise ValueError("La sensibilidad configurada debe incluir solo la curva M1_22A")
    salida.mkdir(parents=True, exist_ok=True)
    total = len(planes) * len(matriz)
    print(f"Matriz: {len(planes)} planes × {len(curvas)} curvas × {len(climas)} climas = {total} corridas")
    print(f"IFC: {ruta_ifc}\nt0 UTC: {t0}; zona local declarada: {dossier['zona_horaria_local']}")

    n = 0
    clima_usado = {}
    for plan in planes:
        vertidas = plan.aplicar(originales)
        grilla = voxelizar_vertidas(vertidas, config["paso_voxel_m"])
        for escenario_id, (curva_cfg, hidratacion), clima_cfg, peso in matriz:
            n += 1
            config_corrida = dict(config)
            overlay = climas_config[str(clima_cfg["id"])]
            for campo in ("clima", "sitio", "entorno", "condicion_superficie"):
                if campo in overlay:
                    config_corrida[campo] = overlay[campo]
            config_corrida.setdefault("clima", {})["inicio_obra_utc"] = t0
            observador = SeguimientoMetricasDecision(
                grilla, float(dossier.get("profundidad_s25_m", 0.025)),
                tuple(dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"])))
            print(f"[{n}/{total}] plan={plan.plan_id}, escenario={escenario_id}")
            resultado_sim, bordes, aire = demo.correr(
                config_corrida, grilla, vertidas, observadores=(observador,),
                modelo_hidratacion=hidratacion)
            if aire is None or aire.es_sintetica:
                raise ValueError(f"Clima {clima_cfg['id']}: se requiere una trayectoria real/reanálisis, no sintética ni constante")
            clima_usado[clima_cfg["id"]] = aire.descripcion
            metricas = observador.resumen()
            resultados.append(ResultadoEscenario(
                plan.plan_id, escenario_id, peso,
                metricas["M_T_c"], metricas["M_delta_c"], metricas["t_M_T_h"],
                metricas["t_M_delta_h"], metricas["xyz_M_T_m"], metricas["xyz_S25_min_m"]))

    # Convergencia en una combinación común de plan/curva/clima; se cambia una resolución
    # a la vez al comparar los resultados, y se conserva el factorial completo para auditoría.
    cfg_conv = dossier.get("convergencia_numerica", {})
    if cfg_conv.get("activo", False):
        plan_ref = next((p for p in planes if p.plan_id == cfg_conv.get("plan_id")), None)
        escenario_ref = next((x for x in matriz if x[0] == cfg_conv.get("escenario_id")), None)
        dxs = [float(x) for x in cfg_conv.get("pasos_voxel_m", [])]
        dts = [float(x) for x in cfg_conv.get("pasos_tiempo_max_s", [])]
        if plan_ref is None or escenario_ref is None or len(set(dxs)) < 2 or len(set(dts)) < 2:
            raise ValueError("convergencia_numerica requiere plan/escenario válidos y al menos dos dx y dos dt")
        escenario_id, (curva_cfg, hidratacion), clima_cfg, _ = escenario_ref
        vertidas_ref = plan_ref.aplicar(originales)
        overlay = climas_config[str(clima_cfg["id"])]
        config_ref = dict(config)
        for campo in ("clima", "sitio", "entorno", "condicion_superficie"):
            if campo in overlay:
                config_ref[campo] = overlay[campo]
        config_ref.setdefault("clima", {})["inicio_obra_utc"] = t0
        filas_conv = []
        for dx in dxs:
            grilla_conv = voxelizar_vertidas(vertidas_ref, dx)
            for dt_max in dts:
                config_corrida = dict(config_ref)
                config_corrida["paso_voxel_m"] = dx
                config_corrida["paso_tiempo_max_s"] = dt_max
                observador = SeguimientoMetricasDecision(
                    grilla_conv, float(dossier.get("profundidad_s25_m", 0.025)),
                    tuple(dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"])))
                print(f"Convergencia: plan={plan_ref.plan_id}, escenario={escenario_id}, dx={dx:g} m, dt_max={dt_max:g} s")
                _, _, aire_conv = demo.correr(config_corrida, grilla_conv, vertidas_ref,
                                               observadores=(observador,), modelo_hidratacion=hidratacion)
                if aire_conv is None or aire_conv.es_sintetica:
                    raise ValueError(f"Clima {clima_cfg['id']}: la convergencia requiere el mismo clima real")
                met = observador.resumen()
                filas_conv.append({"dx_m": dx, "paso_tiempo_max_s": dt_max,
                                   "M_T_c": met["M_T_c"], "t_M_T_h": met["t_M_T_h"],
                                   "M_delta_c": met["M_delta_c"], "t_M_delta_h": met["t_M_delta_h"],
                                   "clasifica_T": met["M_T_c"] > limite_t,
                                   "clasifica_delta": met["M_delta_c"] > limite_delta})
        dx_fino, dt_fino = min(dxs), min(dts)
        ref = next(f for f in filas_conv if f["dx_m"] == dx_fino and f["paso_tiempo_max_s"] == dt_fino)
        for f in filas_conv:
            espacial = next(x for x in filas_conv if x["dx_m"] == dx_fino and
                            x["paso_tiempo_max_s"] == f["paso_tiempo_max_s"])
            temporal = next(x for x in filas_conv if x["dx_m"] == f["dx_m"] and
                            x["paso_tiempo_max_s"] == dt_fino)
            f["delta_M_T_vs_dx_fino_misma_dt_c"] = f["M_T_c"] - espacial["M_T_c"]
            f["delta_M_delta_vs_dx_fino_misma_dt_c"] = f["M_delta_c"] - espacial["M_delta_c"]
            f["delta_M_T_vs_dt_fino_mismo_dx_c"] = f["M_T_c"] - temporal["M_T_c"]
            f["delta_M_delta_vs_dt_fino_mismo_dx_c"] = f["M_delta_c"] - temporal["M_delta_c"]
            f["delta_M_T_vs_mas_fino_c"] = f["M_T_c"] - ref["M_T_c"]
            f["delta_M_delta_vs_mas_fino_c"] = f["M_delta_c"] - ref["M_delta_c"]
        with (salida / "convergencia_numerica.csv").open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.DictWriter(archivo, fieldnames=list(filas_conv[0]))
            escritor.writeheader()
            escritor.writerows(filas_conv)

    # Valida que los pesos usados correspondan a escenarios discretos conjuntos.
    resumen = agregar_probabilidades(resultados, limite_T_c=limite_t, limite_delta_c=limite_delta,
                                     limite_P_f=limite_pf)
    guardar_resultados_escenario_csv(salida / "resultados_por_escenario.csv", resultados,
                                     limite_T_c=limite_t, limite_delta_c=limite_delta)
    with (salida / "probabilidades_por_plan.csv").open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=list(resumen[0]))
        escritor.writeheader()
        escritor.writerows(resumen)
    comparacion = []
    por_plan = {p.plan_id: p for p in planes}
    for fila in resumen:
        comparacion.append({"plan_id": fila["plan_id"], "categoria": categorias[fila["plan_id"]],
            "D_h": por_plan[fila["plan_id"]].duracion_obra_h(
            float(dossier.get("t0_h", 0.0))), "P_T": fila["P_T"], "P_delta": fila["P_delta"],
            "P_f": fila["P_f"], "limite_P_f": limite_pf, "cumple_P_f": fila["cumple_P_f"]})
    factibles = [f for f in comparacion if f["cumple_P_f"]]
    for fila in comparacion:
        fila["mas_rapido_factible_evaluado"] = bool(factibles and fila["plan_id"] == min(factibles, key=lambda x: x["D_h"])["plan_id"])
    with (salida / "comparacion_planes.csv").open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=list(comparacion[0]))
        escritor.writeheader()
        escritor.writerows(comparacion)
    archivos_fuente = [ruta_ifc, ruta_cfg]
    archivos_fuente.extend((RAIZ / c).resolve() for c in dossier.get("config_base", []))
    archivos_fuente.extend((RAIZ / c["archivo_source_terms"]).resolve() for c, _ in curvas)
    archivos_fuente.extend((RAIZ / c["archivo_source_terms"]).resolve()
                            for c in dossier.get("curvas_sensibilidad", []))
    archivos_fuente.extend((RAIZ / c["config_jsonc"]).resolve() for c in climas)
    archivos_fuente.extend(RAIZ / ruta for ruta in (
        "scripts/evaluar_dossier_final.py", "scripts/demo_e10_entorno.py",
        "scripts/preparar_climas_gefs_tres_trayectorias.py",
        "scripts/graficar_climas_gefs_tres_trayectorias.py",
        "scripts/graficar_curvas_hidratacion.py",
        "scripts/diagnostico_horizonte_final.py",
        "scripts/verificar_convergencia_final.py",
        "scripts/demo_pipeline_ifc_3d.py", "src/mcplan/thermal/decision.py",
        "src/mcplan/thermal/hidratacion.py", "src/mcplan/thermal/vertidas.py",
        "src/mcplan/thermal/diffusion_nd.py", "src/mcplan/thermal/bordes.py",
        "src/mcplan/geometry/ifc.py", "src/mcplan/geometry/voxelizacion.py"))
    for c in climas:
        clima_cfg = climas_config[str(c["id"])].get("clima", {})
        snapshot = clima_cfg.get("snapshot", {})
        carpeta = (RAIZ / snapshot["carpeta"]).resolve()
        if not carpeta.is_dir():
            raise FileNotFoundError(f"No existe carpeta local del snapshot climático: {carpeta}")
        archivos_fuente.extend(f for f in carpeta.rglob("*") if f.is_file())
    carpeta_ensamble = RAIZ / "data/weather/2026-11-16T213220Z_GEFS_ensemble"
    archivos_fuente.extend(carpeta_ensamble / nombre for nombre in
                           ("original.json", "ensemble_metadata.json"))
    seleccion_clima = RAIZ / "data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias/selection_metadata.json"
    if seleccion_clima.is_file():
        archivos_fuente.append(seleccion_clima)
    archivos_fuente = sorted(set(archivos_fuente))
    git_metadata_error = None
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=RAIZ, check=True,
                                capture_output=True, text=True).stdout.strip()
        git_status = subprocess.run(["git", "status", "--porcelain"], cwd=RAIZ, check=True,
                                    capture_output=True, text=True).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError) as exc:
        # No informar falsamente que el árbol de trabajo está limpio cuando Git
        # no permite consultar su estado; el commit tampoco identifica el código
        # evaluado si no se puede comprobar esa relación.
        commit = None
        git_status = None
        git_metadata_error = str(exc)
    manifest = {
        "generado_utc": datetime.now(timezone.utc).isoformat(), "commit": commit,
        "tipo_evaluacion": "sensibilidad_curva_M1_22A" if sensibilidad else "nominal_M1_SS",
        "working_tree_status": git_status,
        "git_metadata_error": git_metadata_error,
        "ifc": str(ruta_ifc.relative_to(RAIZ)), "inicio_obra_utc": t0,
        "inicio_obra_local": dossier.get("inicio_obra_local"),
        "zona_horaria_local": dossier["zona_horaria_local"],
        "sitio": config.get("sitio"),
        "curvas_nominales": [c["id"] for c, _ in curvas],
        "curvas_sensibilidad": [c["id"] for c in dossier.get("curvas_sensibilidad", [])],
        "pesos_clima": {str(c["id"]): float(c["peso"]) for c in climas},
        "horizonte_h": config["tiempo_final_h"], "paso_voxel_m": config["paso_voxel_m"],
        "paso_tiempo_max_s": config["paso_tiempo_max_s"],
        "criterios": {"M_T_c": limite_t, "M_delta_c": limite_delta,
                      "probabilidad_fallo_maxima": limite_pf},
        "clima_usado": clima_usado,
        "archivos": {str(p.relative_to(RAIZ)): sha256(p) for p in archivos_fuente},
        "definicion_S25": {"profundidad_m": dossier.get("profundidad_s25_m", 0.025),
                           "caras": dossier.get("caras_s25", ["x-", "x+", "y-", "y+", "z-", "z+"]),
                           "reconstruccion": "temperatura de cara con Robin/Dirichlet/Neumann y resistencia de protecciones; interpolación lineal al interior"},
        "supuesto_activacion": "cada vertida se activa instantáneamente en activaciones_h; las duraciones planificadas calculan D y capacidad, no colocación parcial",
        "supuesto_probabilidad": "pesos de escenarios conjuntos" if dossier.get("escenarios_conjuntos") else "independencia entre curva y clima declarada",
    }
    (salida / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    elegido = min(factibles, key=lambda x: x["D_h"]) if factibles else None
    print(f"\nSalidas: {salida}\nPlan factible más rápido entre los evaluados: "
          f"{elegido['plan_id'] if elegido else 'ninguno'}")
    if not factibles:
        print("No se identificó un plan que cumpla P_f ≤ 0.05 entre los escenarios evaluados.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR DE CONFIGURACIÓN/DATOS: {exc}", file=sys.stderr)
        raise SystemExit(2)
