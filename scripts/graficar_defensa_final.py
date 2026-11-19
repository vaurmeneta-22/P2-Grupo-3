"""Genera figuras listas para la defensa y el HTML desde la corrida congelada.

Lee los resultados ya calculados de ``evaluar_dossier_final.py`` y ejecuta una
serie temporal representativa por plan (M1_SS + clima central) para producir
curvas y cortes 3D. No modifica el solver ni vuelve a calcular la matriz de 24
escenarios.
"""
from __future__ import annotations

import csv
import copy
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(RAIZ / "src"))

import demo_e10_entorno as demo
import demo_pipeline_ifc_3d as pipe
import evaluar_dossier_final as evaluador
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import SeguimientoMetricasDecision
from mcplan.weather.loader import construir_aire


DOSSIER = RAIZ / "configs/defensa_final_planes_100m3h.json"
RESULTADOS = RAIZ / "results/defensa_final_100m3h"
FIGURAS = RAIZ / "reports/Semana Final/evidencia_semana_final/comparacion_100m3h"
ORDEN = [
    "caso_ifc_original_sin_manta", "control_10h_T10_sin_manta",
    "vertida_unica_100m3h", "propuesto_10h_60mm_T10",
    "secuencia_continua_100m3h", "secuencia_continua_12h",
]
ETIQUETAS = {
    "caso_ifc_original_sin_manta": "IFC original\nsin manta",
    "control_10h_T10_sin_manta": "Control 10 h\nsin manta",
    "vertida_unica_100m3h": "Vertida única\n100 m³/h",
    "propuesto_10h_60mm_T10": "Propuesto\n10 h · PUR 60 mm",
    "secuencia_continua_100m3h": "Secuencia rápida\n100 m³/h",
    "secuencia_continua_12h": "Secuencia lenta\n12 h · PUR 60 mm",
}
COLORES = {
    "caso_ifc_original_sin_manta": "#64748b",
    "control_10h_T10_sin_manta": "#94a3b8",
    "vertida_unica_100m3h": "#f59e0b",
    "propuesto_10h_60mm_T10": "#2563eb",
    "secuencia_continua_100m3h": "#0f9d78",
    "secuencia_continua_12h": "#8b5cf6",
}


def leer_csv(ruta: Path) -> list[dict[str, str]]:
    with ruta.open(encoding="utf-8-sig", newline="") as archivo:
        return list(csv.DictReader(archivo))


def preparar_estilo() -> None:
    plt.rcParams.update({
        "figure.facecolor": "#f8fafc", "axes.facecolor": "#ffffff",
        "axes.edgecolor": "#cbd5e1", "axes.labelcolor": "#334155",
        "text.color": "#0f172a", "xtick.color": "#475569",
        "ytick.color": "#475569", "grid.color": "#e2e8f0",
        "grid.alpha": 0.8, "font.family": "DejaVu Sans",
        "axes.titleweight": "bold", "axes.titlesize": 12,
        "axes.labelsize": 10, "legend.frameon": False,
    })


def grafico_resumen(comparacion: list[dict], escenarios: list[dict]) -> None:
    by_plan = {r["plan_id"]: r for r in comparacion}
    grupos: dict[str, list[dict]] = {p: [] for p in ORDEN}
    for fila in escenarios:
        grupos[fila["plan_id"]].append(fila)
    tmax = {p: max(float(x["M_T_c"]) for x in grupos[p]) for p in ORDEN}
    delta = {p: max(float(x["M_delta_c"]) for x in grupos[p]) for p in ORDEN}
    labels = [ETIQUETAS[p].replace("\n", " ") for p in ORDEN]
    colors = [COLORES[p] for p in ORDEN]
    series = [
        ([float(by_plan[p]["D_h"]) for p in ORDEN], None, "Duración de obra", "D [h]"),
        ([100 * float(by_plan[p]["P_f"]) for p in ORDEN], 5, "Probabilidad de fallo", "Pf [%]"),
        ([delta[p] for p in ORDEN], 20, "Diferencial máximo", "ΔT [°C]"),
        ([tmax[p] for p in ORDEN], 70, "Temperatura máxima", "Tmax [°C]"),
    ]
    fig, ejes = plt.subplots(1, 4, figsize=(18, 7.0), sharey=True)
    ypos = np.arange(len(ORDEN))
    for ax, (values, limit, title, xlabel) in zip(ejes, series):
        bars = ax.barh(ypos, values, color=colors, height=0.66)
        if limit is not None:
            ax.axvline(limit, color="#dc2626", linestyle="--", linewidth=1.5,
                       label=f"Límite {limit} {xlabel.split()[-1]}")
            ax.legend(loc="lower right", fontsize=8)
        ax.set_title(title, weight="bold")
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="x")
        ax.set_axisbelow(True)
        for bar, value in zip(bars, values):
            ax.text(bar.get_width() + max(values) * 0.025, bar.get_y() + bar.get_height()/2,
                    f"{value:.2f}" if value < 10 else f"{value:.1f}",
                    va="center", ha="left", fontsize=8, color="#334155")
        ax.set_xlim(0, max(values) * 1.25 if max(values) else 1)
    ejes[0].set_yticks(ypos, labels, fontsize=8)
    ejes[0].invert_yaxis()
    fig.suptitle("Comparación de planes · 24 escenarios ponderados por plan", fontsize=16, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIGURAS / "resumen_planes.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def grafico_matriz(escenarios: list[dict]) -> None:
    by_id = {r["plan_id"]: r for r in escenarios}
    curves = list(dict.fromkeys(r["escenario_id"].split("__")[0] for r in escenarios))
    climates = ["frio", "central", "calido"]
    scenario_ids = [f"{curve}__GEFS_20261116_{climate}" for curve in curves for climate in climates]
    matriz = np.array([[float(next(r for r in escenarios if r["plan_id"] == p and r["escenario_id"] == s)["M_delta_c"])
                        for s in scenario_ids] for p in ORDEN])
    fig, ax = plt.subplots(figsize=(17, 6.2))
    im = ax.imshow(matriz, cmap="YlOrRd", vmin=0, vmax=50, aspect="auto")
    ax.set_yticks(np.arange(len(ORDEN)), [ETIQUETAS[p].replace("\n", " ") for p in ORDEN])
    short = {"frio": "fría", "central": "central", "calido": "cálida"}
    ticklabels = [f"{curve.replace('_SS','').replace('_22A','*')}\n{short[climate]}"
                  for curve in curves for climate in climates]
    ax.set_xticks(np.arange(len(scenario_ids)), ticklabels, rotation=0, ha="center", fontsize=7.3)
    ax.set_title("Máximo ΔT por plan y escenario", pad=12, fontsize=15, weight="bold")
    ax.set_xlabel("Curva de hidratación × trayectoria climática  (* = campaña 22A; sin * = SS)")
    ax.set_ylabel("Plan")
    for i in range(len(ORDEN)):
        for j in range(len(scenario_ids)):
            val = matriz[i, j]
            if val > 20:
                ax.scatter(j, i, marker="x", s=33, linewidth=1.5, color="#7f1d1d")
    cbar = fig.colorbar(im, ax=ax, pad=0.015)
    cbar.set_label("ΔT máximo [°C] · × marca excedencia de 20 °C")
    fig.tight_layout()
    fig.savefig(FIGURAS / "matriz_diferencial_escenarios.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def grafico_cronogramas(dossier: dict) -> None:
    planes = {p["id"]: p for p in dossier["planes"]}
    fig, ax = plt.subplots(figsize=(13, 6.3))
    y = np.arange(len(ORDEN))[::-1]
    paleta_vertidas = {1:"#1d4ed8", 2:"#0f766e", 3:"#d97706", 4:"#7c3aed", 5:"#db2777"}
    for yi, plan_id in zip(y, ORDEN):
        p = planes[plan_id]
        if plan_id == "vertida_unica_100m3h":
            d = float(p["duraciones_h"]["1"])
            ax.barh(yi, d, left=0, height=0.46, color=COLORES[plan_id], edgecolor="white")
            ax.text(d/2, yi, "Todo activado en t=0", color="white", ha="center", va="center", fontsize=8)
        else:
            for i in range(1, 6):
                start = float(p["activaciones_h"][str(i)])
                duration = float(p["duraciones_h"][str(i)])
                ax.barh(yi, duration, left=start, height=0.48, color=paleta_vertidas[i], edgecolor="white", linewidth=0.7)
                if duration >= 0.55:
                    ax.text(start + duration/2, yi, f"V{i}", color="white", ha="center", va="center", fontsize=8, weight="bold")
        fin = max(float(p["activaciones_h"][str(i)]) + float(p["duraciones_h"][str(i)]) for i in range(1, 6))
        ax.text(fin + 0.18, yi, f"D={fin:.2f} h", va="center", fontsize=8, color="#334155")
        manta = next((x["inicio_h"] for grupo in p["entorno"]["protecciones_por_vertida"].values()
                      for x in grupo if x["nombre"] == "manta_PUR"), None)
        if manta is not None:
            ax.scatter([manta], [yi], marker="D", s=35, color="#0f172a", zorder=5)
    ax.set_yticks(y, [ETIQUETAS[p].replace("\n", " ") for p in ORDEN])
    ax.set_xlabel("Horas desde el inicio común · 19 nov 2026, 10:00 Ovalle")
    ax.set_title("Cronogramas comparados · diamante = instalación de manta", fontsize=14, weight="bold")
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(30, max(float(r["D_h"]) for r in leer_csv(RESULTADOS / "comparacion_planes.csv")) + 4))
    leyenda = [plt.Line2D([0],[0], color=c, lw=7, label=f"Vertida {i}") for i,c in paleta_vertidas.items()]
    leyenda.append(plt.Line2D([0],[0], marker="D", linestyle="", color="#0f172a", label="Instala manta"))
    ax.legend(handles=leyenda, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    fig.savefig(FIGURAS / "cronogramas_planes.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def correr_historiales(dossier: dict) -> tuple[dict[str, list[dict]], dict]:
    config, _ = evaluador.cargar_configuracion_dossier(DOSSIER)
    curvas = evaluador.leer_curvas(dossier)
    hidratacion = next(modelo for cfg, modelo in curvas if cfg["id"] == "M1_SS")
    climas = dossier["climas"]
    clima_cfg = next(c for c in climas if c["id"].endswith("_central"))
    clima_overlay = evaluador.leer_jsonc((RAIZ / clima_cfg["config_jsonc"]).resolve())
    clima_local = dict(clima_overlay.get("clima") or {})
    clima_local["inicio_obra_utc"] = dossier["inicio_obra_utc"]
    ruta_ifc = (RAIZ / config["ifc"]).resolve()
    originales = extraer_vertidas(ruta_ifc)
    plans_cfg = {str(p["id"]): p for p in dossier["planes"]}
    history = {}
    salida_campos = {}
    ids_run = [
        "caso_ifc_original_sin_manta", "control_10h_T10_sin_manta",
        "vertida_unica_100m3h", "propuesto_10h_60mm_T10",
        "secuencia_continua_100m3h", "secuencia_continua_12h",
    ]
    config_clima = clima_overlay
    for plan_id in ids_run:
        print(f"Serie térmica representativa: {plan_id} · M1_SS · clima central")
        pcfg = plans_cfg[plan_id]
        plan = evaluador.cargar_plan(pcfg, originales)
        vertidas = plan.aplicar(originales)
        grilla = voxelizar_vertidas(vertidas, float(dossier["paso_voxel_m"]))
        cfg_run = dict(config)
        for key in ("clima", "sitio", "entorno", "condicion_superficie"):
            if key in config_clima:
                cfg_run[key] = config_clima[key]
        cfg_run["clima"] = dict(clima_local)
        cfg_run["entorno"] = copy.deepcopy(pcfg["entorno"])
        observador = SeguimientoMetricasDecision(
            grilla, float(dossier["profundidad_s25_m"]), tuple(dossier["caras_s25"]),
            guardar_historial=True)
        resultado, _, aire = demo.correr(
            cfg_run, grilla, vertidas, observadores=(observador,), modelo_hidratacion=hidratacion)
        if aire is None or aire.es_sintetica:
            raise ValueError("La serie representativa requiere el snapshot climático real")
        history[plan_id] = observador.historial
        if plan_id == "propuesto_10h_60mm_T10":
            salida_campos = {"grilla": grilla, "resultado": resultado, "resumen": observador.resumen()}
    ruta_csv = FIGURAS / "series_temporales_representativas.csv"
    with ruta_csv.open("w", newline="", encoding="utf-8") as archivo:
        campos = ["plan_id", "tiempo_h", "temperatura_maxima_c", "temperatura_media_c",
                  "temperatura_minima_voxel_c", "temperatura_minima_S25_c", "diferencia_maxima_S25_c"]
        writer = csv.DictWriter(archivo, fieldnames=campos)
        writer.writeheader()
        for plan_id, items in history.items():
            for item in items:
                writer.writerow({"plan_id": plan_id, **item})
    return history, salida_campos


def grafico_temperaturas(history: dict[str, list[dict]]) -> None:
    fig, ejes = plt.subplots(1, 3, figsize=(18, 5.2))
    for plan_id in ORDEN:
        h = history[plan_id]
        t = [x["tiempo_h"] for x in h]
        color = COLORES[plan_id]
        ejes[0].plot(t, [x["temperatura_maxima_c"] for x in h], color=color, lw=1.7,
                     label=ETIQUETAS[plan_id].replace("\n", " "))
        ejes[1].plot(t, [x["temperatura_media_c"] for x in h], color=color, lw=1.6,
                     label=ETIQUETAS[plan_id].replace("\n", " "))
        ejes[1].plot(t, [x["temperatura_minima_voxel_c"] for x in h], color=color, lw=0.9, alpha=0.48)
        ejes[2].plot(t, [x["diferencia_maxima_S25_c"] for x in h], color=color, lw=1.7,
                     label=ETIQUETAS[plan_id].replace("\n", " "))
    ejes[0].axhline(70, color="#dc2626", ls="--", lw=1.3)
    ejes[0].set_title("Máxima global")
    ejes[0].set_ylabel("Temperatura [°C]")
    ejes[1].set_title("Media (línea) y mínimo voxel (tenue)")
    ejes[1].set_ylabel("Temperatura [°C]")
    ejes[2].axhline(20, color="#dc2626", ls="--", lw=1.3, label="Límite 20 °C")
    ejes[2].set_title("ΔT = Tmax global − Tmin S25")
    ejes[2].set_ylabel("Diferencia [°C]")
    for ax in ejes:
        ax.set_xlabel("Tiempo desde inicio de obra [h]")
        ax.grid(True)
    ejes[0].legend(fontsize=7, loc="best")
    ejes[2].legend(fontsize=7, loc="best")
    fig.suptitle("Evolución térmica · M1_SS + GEFS central", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGURAS / "evolucion_temperaturas_y_delta.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def grafico_cortes(campos: dict) -> None:
    grilla, resultado = campos["grilla"], campos["resultado"]
    xyz = campos["resumen"]["xyz_M_T_m"]
    iz = int(np.clip(round((xyz[2] - grilla.origen_m[2]) / grilla.paso_m[2] - 0.5), 0, grilla.forma[2]-1))
    tiempos = [24, 72, 168, 336]
    indices = [int(np.argmin(np.abs(resultado.tiempos_s / 3600.0 - t))) for t in tiempos]
    x0, y0 = grilla.origen_m[:2]
    dx, dy = grilla.paso_m[:2]
    extent = (x0, x0 + grilla.forma[0]*dx, y0, y0 + grilla.forma[1]*dy)
    fig, ejes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    imagen = None
    for ax, t, idx in zip(ejes.flat, tiempos, indices):
        field = resultado.temperaturas_c[idx][:, :, iz]
        active = resultado.mascaras_activas[idx][:, :, iz]
        view = np.ma.masked_where(~active, field).T
        imagen = ax.imshow(view, origin="lower", extent=extent, cmap="inferno", vmin=10, vmax=70,
                           interpolation="nearest", aspect="equal")
        ax.set_title(f"t = {resultado.tiempos_s[idx]/3600:.0f} h")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    fig.colorbar(imagen, ax=ejes.ravel().tolist(), label="Temperatura [°C]", shrink=0.82)
    fig.suptitle(f"Corte horizontal · plan 60 mm · capa z≈{grilla.origen_m[2] + (iz+0.5)*grilla.paso_m[2]:.2f} m",
                 fontsize=14, weight="bold")
    fig.savefig(FIGURAS / "cortes_temperatura_60mm.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def grafico_entradas(dossier: dict) -> None:
    # Ocho curvas de fuente del archivo original del curso.
    ruta_fuente = RAIZ / dossier["curvas"][0]["archivo_source_terms"]
    lineas = ruta_fuente.read_text(encoding="utf-8").splitlines()
    encabezado = next(line.lstrip()[1:].strip() for line in lineas
                      if line.lstrip().startswith("#") and "time_h" in line and "M1_SS" in line)
    nombres = encabezado.split()
    datos = [line for line in lineas if line.strip() and not line.lstrip().startswith("#")]
    valores = np.genfromtxt(datos, dtype=float)
    fig, ejes = plt.subplots(1, 2, figsize=(14, 5.3))
    palette = plt.cm.viridis(np.linspace(0.08, 0.92, len(nombres)-1))
    for c, name in zip(palette, nombres[1:]):
        ejes[0].plot(valores[:, 0], valores[:, nombres.index(name)], label=name, lw=1.7, color=c)
    ejes[0].set_title("Fuente volumétrica de hidratación")
    ejes[0].set_xlabel("Edad de la mezcla [h]")
    ejes[0].set_ylabel("Q̇ [W/m³]")
    ejes[0].grid(True)
    ejes[0].legend(ncol=2, fontsize=8)

    t0 = dossier["inicio_obra_utc"]
    for clima in dossier["climas"]:
        overlay = evaluador.leer_jsonc((RAIZ / clima["config_jsonc"]).resolve())
        cfg = dict(overlay.get("clima") or {})
        cfg["inicio_obra_utc"] = t0
        aire = construir_aire(cfg, tiempo_final_s=336*3600, raiz_repo=RAIZ)
        serie = aire.serie
        x = serie.segundos_desde(__import__("mcplan.weather.serie", fromlist=["a_utc"]).a_utc(t0)) / 3600
        ax = ejes[1]
        label = "fría" if "frio" in clima["id"] else ("cálida" if "calido" in clima["id"] else "central")
        ax.plot(x, serie.temperatura_c, label=label, lw=1.7)
    ejes[1].set_title("Temperatura ambiente · clima GEFS")
    ejes[1].set_xlabel("Tiempo desde inicio de obra [h]")
    ejes[1].set_ylabel("Aire [°C]")
    ejes[1].grid(True)
    ejes[1].legend()
    fig.suptitle("Datos de entrada", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGURAS / "curvas_hidratacion_y_clima.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    # Exporta además cada panel por separado para la presentación por pantallas.
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    palette = plt.cm.viridis(np.linspace(0.06, 0.94, len(nombres) - 1))
    for color, name in zip(palette, nombres[1:]):
        ax.plot(valores[:, 0], valores[:, nombres.index(name)], label=name,
                lw=2.5 if name == "M1_SS" else 1.8, color=color)
    ax.set_title("Fuente volumétrica de hidratación", weight="bold")
    ax.set_xlabel("Edad desde la colocación [h]")
    ax.set_ylabel("q̇ [W/m³]")
    ax.grid(True)
    ax.legend(ncol=2, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURAS / "fuente_hidratacion_8_curvas.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for clima in dossier["climas"]:
        overlay = evaluador.leer_jsonc((RAIZ / clima["config_jsonc"]).resolve())
        cfg = dict(overlay.get("clima") or {})
        cfg["inicio_obra_utc"] = dossier["inicio_obra_utc"]
        aire = construir_aire(cfg, tiempo_final_s=336*3600, raiz_repo=RAIZ)
        serie = aire.serie
        x = serie.segundos_desde(__import__("mcplan.weather.serie", fromlist=["a_utc"]).a_utc(dossier["inicio_obra_utc"])) / 3600
        label = "Fría" if "frio" in clima["id"] else ("Cálida" if "calido" in clima["id"] else "Central")
        ax.plot(x, serie.temperatura_c, label=label, lw=2.0)
    ax.set_title("Temperatura ambiente · clima GEFS", weight="bold")
    ax.set_xlabel("Tiempo desde el inicio de obra [h]")
    ax.set_ylabel("Temperatura del aire [°C]")
    ax.set_xlim(0, 336)
    ax.grid(True)
    ax.legend(ncol=3, loc="best")
    fig.tight_layout()
    fig.savefig(FIGURAS / "temperatura_aire_gefs.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    FIGURAS.mkdir(parents=True, exist_ok=True)
    preparar_estilo()
    dossier = evaluador.leer_jsonc(DOSSIER)
    comparacion = leer_csv(RESULTADOS / "comparacion_planes.csv")
    escenarios = leer_csv(RESULTADOS / "resultados_por_escenario.csv")
    if len(escenarios) != len(ORDEN) * 24:
        raise ValueError(f"Se esperaban {len(ORDEN)*24} resultados y hay {len(escenarios)}")
    grafico_resumen(comparacion, escenarios)
    grafico_matriz(escenarios)
    grafico_cronogramas(dossier)
    grafico_entradas(dossier)

    config, _ = evaluador.cargar_configuracion_dossier(DOSSIER)
    vertidas = extraer_vertidas((RAIZ / config["ifc"]).resolve())
    grilla = voxelizar_vertidas(vertidas, float(dossier["paso_voxel_m"]))
    pipe.dibujar_geometria_ifc_ensamblada(vertidas, FIGURAS / "geometria_fundacion_ifc.png")
    pipe.dibujar_geometria(grilla, None, FIGURAS / "geometria_voxelizada_por_vertida.png")

    history, campos = correr_historiales(dossier)
    grafico_temperaturas(history)
    grafico_cortes(campos)
    print(f"Figuras y series guardadas en: {FIGURAS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
