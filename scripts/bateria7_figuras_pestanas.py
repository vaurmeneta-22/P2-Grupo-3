"""Figuras para las pestañas «Caso base», «Plan 1» y «Alternativas» de la presentación.

Solo grafica resultados ya calculados (results/bateria7/cache y campos); no ejecuta simulaciones.
Las escalas de color 3D se ajustan al rango real de los campos de cada grupo de comparación
(mínimo y máximo de temperatura de todas las vistas del grupo, redondeados a °C enteros) y se
comparten dentro del grupo. Escribe además reports/Defensa/bateria7/tablas/escalas_3d.json.
"""
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
OUT = RAIZ / "reports/Defensa/bateria7/figuras"
CFG = json.loads((RAIZ / "configs/defensa/bateria7_casos.json").read_text(encoding="utf-8"))
PLANES = {p["id"]: p for p in CFG["planes"]}
FRIO, CENTRAL = "M4_SS__GEFS_20261116_frio", "M1_SS__GEFS_20261116_central"
PCOL = ["#7be0bd", "#79c8bb", "#ffc17f", "#a9a0ff", "#f38ba8"]
plt.rcParams.update({"font.size": 10})


def serie(plan, e=FRIO):
    return json.loads((SAL / "cache" / f"{plan}__{e}.json").read_text(encoding="utf-8"))["serie"]


def campo(plan, e=FRIO):
    return np.load(SAL / "campos" / f"{plan}__{e}.npz")


def temps_en(plan, t_obj, e=FRIO):
    d = campo(plan, e)
    i = int(np.argmin(np.abs(d["t_h"] - t_obj)))
    return d, i, d["T"][i]


def rango(items, tiempos):
    lo, hi = [], []
    for plan in items:
        for t in tiempos:
            _, _, T = temps_en(plan, t)
            lo.append(float(np.nanmin(T))); hi.append(float(np.nanmax(T)))
    return math.floor(min(lo)), math.ceil(max(hi))


def vista3d(plan, tiempos, titulo, nombre, vmin, vmax):
    d = campo(plan); idx = np.argwhere(d["mascara"]); xyz = d["origen"] + (idx + 0.5) * d["paso"]
    sel = xyz[:, 1] <= 0.5 * (xyz[:, 1].min() + xyz[:, 1].max())
    fig = plt.figure(figsize=(6.5 * len(tiempos), 4.9))
    for j, tobj in enumerate(tiempos):
        i = int(np.argmin(np.abs(d["t_h"] - tobj))); T = d["T"][i]
        ax = fig.add_subplot(1, len(tiempos), j + 1, projection="3d")
        sc = ax.scatter(xyz[sel, 0], xyz[sel, 1], xyz[sel, 2], c=T[sel], cmap="inferno", vmin=vmin, vmax=vmax, s=6, marker="s", linewidths=0)
        ax.set_title(f"t = {d['t_h'][i]:.0f} h · Tmax {np.nanmax(T):.1f} °C", fontsize=10); ax.view_init(24, 125)
        ax.set_box_aspect((np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]) / 2 + 1e-3, np.ptp(xyz[:, 2])))
        ax.set_xticks([-8, -4, 0, 4, 8]); ax.set_yticks([-8, -4, 0]); ax.set_zticks([0, 1, 2, 3, 4]); ax.tick_params(labelsize=8, pad=-2)
        ax.set_xlabel("x [m]", fontsize=9, labelpad=-4); ax.set_ylabel("y [m]", fontsize=9, labelpad=-4); ax.set_zlabel("z [m]", fontsize=9, labelpad=-4)
    cax = fig.add_axes([0.3, 0.13, 0.4, 0.03]); fig.colorbar(sc, cax=cax, orientation="horizontal", label=f"Temperatura [°C] (escala común {vmin}–{vmax} °C)")
    fig.suptitle(f"{titulo} — corte longitudinal, escenario M4_SS / clima frío", fontsize=11); fig.savefig(OUT / nombre, dpi=120); plt.close(fig)


def curvas3(plan, titulo, nombre, col):
    s = serie(plan); t = s["tiempo_h"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.1), sharex=True)
    ax[0].plot(t, s["temperatura_maxima_c"], color=col, lw=2); ax[0].axhline(70, color="k"); ax[0].set_title("Tmax [°C] (límite 70)"); ax[0].set_ylim(0, 75)
    ax[1].plot(t, s["temperatura_minima_S25_c"], color=col, lw=2); ax[1].set_title("Tmin S25 [°C]")
    ax[2].plot(t, s["diferencia_maxima_S25_c"], color=col, lw=2); ax[2].axhline(20, color="k"); ax[2].set_title("ΔT = Tmax − Tmin S25 [°C] (límite 20)")
    y = np.asarray(s["diferencia_maxima_S25_c"]); im = int(np.argmax(y)); ax[2].plot(t[im], y[im], "o", color="k")
    ax[2].annotate(f"máx {y[im]:.2f} °C a {t[im]:.0f} h", (t[im], y[im]), xytext=(-140, -30), textcoords="offset points", fontsize=9)
    for a in ax: a.set_xlabel("Tiempo desde t0 [h]"); a.grid(alpha=.3); a.set_xlim(0, 336)
    ax[0].set_ylabel("Temperatura [°C]"); ax[1].set_ylabel("Temperatura [°C]"); ax[2].set_ylabel("Diferencia térmica ΔT [°C]")
    fig.suptitle(f"{titulo} — escenario M4_SS / clima frío (el de mayor ΔT)"); fig.tight_layout(); fig.savefig(OUT / nombre, dpi=130); plt.close(fig)


def superpos3(items, titulo, nombre, e, solo_delta=False):
    n = 1 if solo_delta else 3
    fig, ax = plt.subplots(1, n, figsize=(9.5 if solo_delta else 16, 4.6 if solo_delta else 4.5)); ax = np.atleast_1d(ax)
    for plan, lab, col, ls, lw in items:
        s = serie(plan, e); t = s["tiempo_h"]
        if solo_delta:
            ax[0].plot(t, s["diferencia_maxima_S25_c"], color=col, ls=ls, lw=lw, label=lab)
        else:
            ax[0].plot(t, s["temperatura_maxima_c"], color=col, ls=ls, lw=lw, label=lab); ax[1].plot(t, s["temperatura_minima_S25_c"], color=col, ls=ls, lw=lw); ax[2].plot(t, s["diferencia_maxima_S25_c"], color=col, ls=ls, lw=lw)
    d = ax[-1]; d.axhline(20, color="k", lw=1.5); d.text(2, 20.5, "límite ΔT = 20 °C", fontsize=9)
    if not solo_delta:
        ax[0].axhline(70, color="k"); ax[0].set_title("Tmax [°C] (límite 70)"); ax[0].set_ylim(0, 75); ax[1].set_title("Tmin S25 [°C]")
    d.set_title("ΔT = Tmax − Tmin S25 [°C] (límite 20)")
    for a in ax: a.set_xlabel("Tiempo desde t0 [h]"); a.grid(alpha=.3); a.set_xlim(0, 336)
    if solo_delta: ax[0].set_ylabel("Diferencia térmica ΔT [°C]")
    else: ax[0].set_ylabel("Temperatura [°C]"); ax[1].set_ylabel("Temperatura [°C]"); ax[2].set_ylabel("Diferencia térmica ΔT [°C]")
    ax[0].legend(fontsize=8, loc="lower right"); fig.suptitle(f"{titulo} — escenario {e.replace('__GEFS_20261116_', ' / clima ')}"); fig.tight_layout(); fig.savefig(OUT / nombre, dpi=130); plt.close(fig)


escalas = {}

# ---------- Caso base (se regenera solo el 3D con escala real) ----------
BASE = [("v1", "bat1_base_original", "1 · Original IFC, 18 °C"), ("v2", "bat2_base_simultaneo", "2 · Simultáneo (t=0), 18 °C"),
        ("v3", "bat3_base_largo_24h", "3 · Largo 24 h, 18 °C"), ("v4", "bat11_base_original_T10", "4 · Original IFC a 10 °C")]
tb = (162, 336)
vmin, vmax = rango([b[1] for b in BASE], tb); escalas["base"] = {"vmin": vmin, "vmax": vmax, "tiempos_h": list(tb)}
for k, plan, nombre in BASE:
    vista3d(plan, tb, nombre, f"base_{k}_3d.png", vmin, vmax)

# ---------- Plan 1: temperatura de colocación ----------
TEMPS = [(10, "bat4_plan1_original_10h", "#1f77b4"), (11, "bat4_T11", "#2ca02c"), (12, "bat4_T12", "#ff7f0e"), (13, "bat4_T13", "#d62728")]
tp = (216, 336)
vmin, vmax = rango([t[1] for t in TEMPS], tp); escalas["plan1_temperaturas"] = {"vmin": vmin, "vmax": vmax, "tiempos_h": list(tp)}
for T, plan, col in TEMPS:
    curvas3(plan, f"Plan 1 con colocación a {T} °C", f"plan1_T{T}_curvas.png", col)
    vista3d(plan, tp, f"Plan 1 · colocación a {T} °C", f"plan1_T{T}_3d.png", vmin, vmax)
for etiq, e in (("frio", FRIO), ("central", CENTRAL)):
    superpos3([(plan, f"{T} °C", col, "-", 2) for T, plan, col in TEMPS], "Plan 1 · temperaturas de colocación 10–13 °C", f"plan1_superposicion_{etiq}.png", e)

# ---------- Alternativas: tres planes con PUR 60 mm a 10 °C ----------
ALT = [("a1", "bat5_plan_simultaneo", "Activación simultánea idealizada (t=0)", "#2ca02c"), ("a2", "bat4_plan1_original_10h", "Programa original de 10 h", "#1f77b4"),
       ("a3", "bat6_plan_largo_24h", "Colocación larga de 24 h", "#d62728")]
vmin, vmax = rango([a[1] for a in ALT] + ["bat7_control_sin_manta"], tp); escalas["alternativas"] = {"vmin": vmin, "vmax": vmax, "tiempos_h": list(tp)}
manta_inst = {"a1": 11.36776, "a2": 15.0, "a3": 29.0}
for k, plan, nombre, col in ALT:
    curvas3(plan, nombre, f"alt_{k}_curvas.png", col)
    vista3d(plan, tp, nombre + " · PUR 60 mm, 10 °C", f"alt_{k}_3d.png", vmin, vmax)
    p = PLANES[plan]; act = {int(i): v for i, v in p["activaciones_h"].items()}; dur = {int(i): v for i, v in p["duraciones_h"].items()}
    simult = k == "a1"
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for i in range(1, 6):
        ax.barh(6 - i, dur[i], left=act[i], color=PCOL[i - 1], edgecolor="k", linewidth=.7, hatch="//" if simult else None, alpha=.55 if simult else 1)
        ax.plot(act[i], 6 - i, "v", color="k", ms=7)
        ax.text(act[i] + dur[i] + .6, 6 - i, f"V{i}: inicio {act[i]:.2f} h" + ("" if simult else f", dur. {dur[i]:.2f} h"), va="center", fontsize=8)
        if i in (1, 2, 5): ax.plot(act[i] + 72, 6 - i, "s", color="#c0392b", ms=6)
    ax.axvline(manta_inst[k], color="#e67e22", ls="--"); ax.plot(manta_inst[k], 0.35, "D", color="#e67e22", ms=8)
    ax.text(manta_inst[k] + .8, 0.35, f"manta PUR instalada a las {manta_inst[k]:.2f} h".replace(".", ",") if False else f"manta PUR a las {manta_inst[k]:.2f} h", fontsize=8, va="center")
    ax.set_yticks(range(1, 6)); ax.set_yticklabels([f"V{i}" for i in (5, 4, 3, 2, 1)]); ax.set_ylim(0, 5.6); ax.set_xlim(0, 100); ax.set_xlabel("horas desde t0 (19-11-2026 10:00)"); ax.grid(axis="x", alpha=.3)
    sub = "▼ activación térmica instantánea al inicio de cada vertida · ◆ instalación de la manta (D + 5 h) · ■ retiro del moldaje (72 h tras colocar; V1, V2, V5)"
    if simult: sub += "\nbarras rayadas: cota por capacidad (6,368 h) usada solo para fijar la manta; la activación simultánea es una idealización térmica, no una duración constructiva de 0 h"
    ax.set_title(f"Cronograma — {nombre}", fontsize=11); fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.text(0.01, 0.01, sub, fontsize=7.5, va="bottom")
    fig.savefig(OUT / f"alt_{k}_cronograma.png", dpi=130); plt.close(fig)
items = [(a[1], a[2], a[3], "-", 2.2) for a in ALT] + [("bat7_control_sin_manta", "Control: 10 h sin manta", "#9467bd", ":", 2.2)]
for etiq, e in (("frio", FRIO), ("central", CENTRAL)):
    superpos3(items, "Alternativas con PUR 60 mm y control sin manta", f"alt_superposicion_{etiq}.png", e, solo_delta=True)
superpos3([("bat7_control_sin_manta", "Control: 10 h sin manta", "#9467bd", ":", 2.4), ("bat4_plan1_original_10h", "Plan 1: 10 h con PUR 60 mm", "#1f77b4", "-", 2.4)],
          "Efecto exclusivo de la manta (mismo programa, 10 °C y moldaje)", "alt_efecto_manta.png", FRIO)
vista3d("bat7_control_sin_manta", tp, "Control: 10 h sin manta, 10 °C", "alt_control_3d.png", vmin, vmax)

(RAIZ / "reports/Defensa/bateria7/tablas/escalas_3d.json").write_text(json.dumps(escalas, indent=1), encoding="utf-8")
print(escalas)
