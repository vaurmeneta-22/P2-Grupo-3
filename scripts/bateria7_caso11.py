"""Caso 11 (base original IFC a 10 °C) y comparaciones 1-11, 11-7, 7-4."""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parents[1]
SAL = RAIZ / "results/bateria7"
PUB = RAIZ / "reports/Semana Final/evidencia_semana_final/comparacion_100m3h/resultados/resultados_por_escenario.csv"
CURVAS = ["M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"]
WC = {"frio": 0.354838709677, "central": 0.322580645161, "calido": 0.322580645161}
CASOS = [("1 Base original (18 °C)", "caso_ifc_original_sin_manta", "bat1_base_original"), ("11 Base original a 10 °C", None, "bat11_base_original_T10"),
         ("7 Control sin manta (10 h, 10 °C)", "control_10h_T10_sin_manta", "bat7_control_sin_manta"), ("4 Plan 1 (10 h, 10 °C, PUR 60 mm)", "propuesto_10h_60mm_T10", "bat4_plan1_original_10h")]
pub = {(r["plan_id"], r["escenario_id"]): r for r in csv.DictReader(PUB.open(encoding="utf-8"))}
def esc(c, w): return f"{c}__GEFS_20261116_{w}"
D = {}
for n, orig, bat in CASOS:
    d = {}
    for c in CURVAS:
        for w in WC:
            e = esc(c, w)
            if orig:
                r = pub[(orig, e)]
                d[e] = dict(MT=float(r["M_T_c"]), MD=float(r["M_delta_c"]), tT=float(r["t_M_T_h"]), tD=float(r["t_M_delta_h"]),
                            xT=(float(r["x_M_T_m"]), float(r["y_M_T_m"]), float(r["z_M_T_m"])), xS=(float(r["x_S25_min_m"]), float(r["y_S25_min_m"]), float(r["z_S25_min_m"])))
            else:
                m = json.loads((SAL / "cache" / f"{bat}__{e}.json").read_text(encoding="utf-8"))["metricas"]
                d[e] = dict(MT=m["M_T_c"], MD=m["M_delta_c"], tT=m["t_M_T_h"], tD=m["t_M_delta_h"], xT=tuple(m["xyz_M_T_m"]), xS=tuple(m["xyz_S25_min_m"]))
    D[n] = d
def w(e): return WC[e.split("_")[-1]] / 8
filas = []
for n, _, bat in CASOS:
    d = D[n]
    eD = max(d, key=lambda e: d[e]["MD"]); eT = max(d, key=lambda e: d[e]["MT"])
    pT = sum(w(e) for e in d if d[e]["MT"] > 70); pD = sum(w(e) for e in d if d[e]["MD"] > 20)
    pF = sum(w(e) for e in d if d[e]["MT"] > 70 or d[e]["MD"] > 20)
    nf = sum(1 for e in d if d[e]["MD"] > 20)
    nend = sum(1 for e in d if d[e]["tD"] >= 312)
    filas.append(dict(caso=n, M_T_max_c=round(d[eT]["MT"], 3), margen_T_c=round(70 - d[eT]["MT"], 3), esc_M_T=eT, t_M_T_h=round(d[eT]["tT"], 1),
                      xyz_M_T=tuple(round(v, 2) for v in d[eT]["xT"]), M_delta_max_c=round(d[eD]["MD"], 3), margen_delta_c=round(20 - d[eD]["MD"], 3),
                      esc_M_delta=eD, t_M_delta_h=round(d[eD]["tD"], 1), xyz_S25_min=tuple(round(v, 2) for v in d[eD]["xS"]),
                      M_delta_min_entre_escenarios_c=round(min(v["MD"] for v in d.values()), 3), escenarios_fallan_delta=nf, P_T=pT, P_delta=pD, P_f=pF,
                      escenarios_con_pico_ultimas_24h=nend))
with (SAL / "tablas/caso11_comparacion.csv").open("w", newline="", encoding="utf-8") as f:
    wr = csv.DictWriter(f, fieldnames=list(filas[0])); wr.writeheader(); wr.writerows(filas)
for r in filas: print(r)
with (SAL / "tablas/caso11_por_escenario.csv").open("w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(["escenario", "peso", "MT_caso1", "MD_caso1", "MT_caso11", "MD_caso11", "t_MD_caso11_h", "MT_caso7", "MD_caso7", "MD_caso4", "dMD_1a11", "dMD_11a7", "dMD_7a4"])
    for e in D[CASOS[0][0]]:
        a, b, c7, c4 = (D[x[0]][e] for x in CASOS)
        wr.writerow([e, f"{w(e):.6f}", f"{a['MT']:.3f}", f"{a['MD']:.3f}", f"{b['MT']:.3f}", f"{b['MD']:.3f}", f"{b['tD']:.1f}", f"{c7['MT']:.3f}", f"{c7['MD']:.3f}", f"{c4['MD']:.3f}",
                     f"{b['MD']-a['MD']:.3f}", f"{c7['MD']-b['MD']:.3f}", f"{c4['MD']-c7['MD']:.3f}"])
print("efectos sobre el máximo de ΔT/Tmax y rango por escenario:")
n1, n11, n7, n4 = (c[0] for c in CASOS)
for a, b, t in ((n1, n11, "temperatura 18→10 (1→11)"), (n11, n7, "secuencia IFC→10 h (11→7)"), (n7, n4, "manta (7→4)")):
    dm = [D[b][e]["MD"] - D[a][e]["MD"] for e in D[a]]; dt = [D[b][e]["MT"] - D[a][e]["MT"] for e in D[a]]
    print(f"{t}: d(max MD)={max(D[b][e]['MD'] for e in D[b]) - max(D[a][e]['MD'] for e in D[a]):+.3f}  por escenario MD {min(dm):+.2f}..{max(dm):+.2f} (media {np.mean(dm):+.2f}); MT {min(dt):+.2f}..{max(dt):+.2f}")
cols = ["#7f7f7f", "#e377c2", "#9467bd", "#1f77b4"]; ls = ["--", "--", ":", "-"]
for e in ("M4_SS__GEFS_20261116_frio", "M4_SS__GEFS_20261116_central", "M1_SS__GEFS_20261116_central"):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True)
    for (n, _, bat), c, l in zip(CASOS, cols, ls):
        s = json.loads((SAL / "cache" / f"{bat}__{e}.json").read_text(encoding="utf-8"))["serie"]; t = s["tiempo_h"]
        ax[0].plot(t, s["temperatura_maxima_c"], color=c, ls=l, label=n); ax[1].plot(t, s["temperatura_minima_S25_c"], color=c, ls=l); ax[2].plot(t, s["diferencia_maxima_S25_c"], color=c, ls=l)
    ax[0].axhline(70, color="k"); ax[2].axhline(20, color="k"); ax[0].legend(fontsize=7, loc="lower right")
    ax[0].set_title("Tmax [°C] (límite 70)"); ax[1].set_title("Tmin S25 [°C]"); ax[2].set_title("ΔT [°C] (límite 20)")
    for a in ax: a.set_xlabel("h desde t0"); a.grid(alpha=.3); a.set_xlim(0, 336)
    fig.suptitle(f"Casos 1, 11, 7 y 4 — escenario {e}"); fig.tight_layout(); fig.savefig(SAL / "figuras" / f"curvas_caso11_{e}.png", dpi=130); plt.close(fig)


def vista4(e, t_obj, nombre):
    fig = plt.figure(figsize=(16, 5.2))
    for k, (n, _, bat) in enumerate(CASOS):
        d = np.load(SAL / "campos" / f"{bat}__{e}.npz")
        i = int(np.argmin(np.abs(d["t_h"] - t_obj))); idx = np.argwhere(d["mascara"]); T = d["T"][i]
        xyz = d["origen"] + (idx + 0.5) * d["paso"]; sel = xyz[:, 1] <= 0.5 * (xyz[:, 1].min() + xyz[:, 1].max())
        ax = fig.add_subplot(1, 4, k + 1, projection="3d")
        sc = ax.scatter(xyz[sel, 0], xyz[sel, 1], xyz[sel, 2], c=T[sel], cmap="inferno", vmin=10, vmax=60, s=5, marker="s", linewidths=0)
        ax.set_title(f"{n}\nt={d['t_h'][i]:.0f} h, máx {np.nanmax(T):.1f} °C", fontsize=8); ax.view_init(24, 125)
        ax.set_box_aspect((np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]) / 2 + 1e-3, np.ptp(xyz[:, 2]))); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    cax = fig.add_axes([0.3, 0.06, 0.4, 0.03]); fig.colorbar(sc, cax=cax, orientation="horizontal", label="T [°C] (escala común 10–60)")
    fig.suptitle(f"{e}: corte longitudinal, t común"); fig.savefig(SAL / "figuras" / nombre, dpi=120); plt.close(fig)


vista4("M4_SS__GEFS_20261116_frio", 164, "3d_caso11_M4_SS_frio_164h.png")
vista4("M4_SS__GEFS_20261116_frio", 336, "3d_caso11_M4_SS_frio_336h.png")
