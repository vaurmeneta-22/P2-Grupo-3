"""Figuras de la defensa a partir del caché (sin recalcular ni usar internet).

    python scripts/defensa_figuras.py base                 # curvas de hidratación y clima
    python scripts/defensa_figuras.py gantt  <ids,...>      # cronogramas
    python scripts/defensa_figuras.py series <id> <curva> <clima> [<id2> ...]
    python scripts/defensa_figuras.py pareto <tabla.csv>    # D vs P_f
    python scripts/defensa_figuras.py cortes <archivo.npz> <salida.png>
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "scripts"))
CACHE = RAIZ / "results/defensa/cache"
FIG = RAIZ / "results/defensa/figuras"
FIG.mkdir(parents=True, exist_ok=True)


def corrida(plan_id, curva, clima, horizonte=1344.0, dx=0.25, dt=600.0):
    mejores = []
    for f in CACHE.glob(f"{plan_id}__{curva}__{clima}__*.json"):
        d = json.loads(f.read_text(encoding="utf-8"))
        if abs(d["horizonte_h"] - horizonte) < 1e-6 and abs(d["dx_m"] - dx) < 1e-9 and abs(d["dt_max_s"] - dt) < 1e-9:
            mejores.append((f.stat().st_mtime, d))
    if not mejores:
        raise FileNotFoundError(f"sin caché para {plan_id} {curva} {clima} h={horizonte}")
    return max(mejores, key=lambda x: x[0])[1]


def fig_base():
    from mcplan.thermal.hidratacion import ModeloFuenteTabulada
    fuente = RAIZ / "data/hidratacion/source_terms.txt"
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for col in ("M1_SS", "M2_SS", "M3_SS", "M4_SS", "M1_22A", "M2_22A", "M3_22A", "M4_22A"):
        m = ModeloFuenteTabulada.desde_source_terms(fuente, col)
        ls = "-" if col.endswith("SS") else "--"
        ax[0].plot(np.asarray(m.tiempos_s) / 3600, m.qdot_w_m3, ls, label=col, lw=2.2 if col == "M1_SS" else 1.2)
    ax[0].set(xlabel="edad desde la colocación [h]", ylabel="Qdot [W/m³]",
              title="Las ocho curvas Qdot del profesor (source_terms.txt)")
    ax[0].legend(ncol=2, fontsize=8); ax[0].grid(alpha=0.3)
    for et, nombre in (("fria", "fría (11/31)"), ("central", "central (10/31)"), ("calida", "cálida (10/31)")):
        d = json.loads((RAIZ / f"data/weather/2026-11-16T213220Z_GEFS_tres_trayectorias_28d/{et}/forecast.json").read_text())
        ax[1].plot(np.arange(len(d["temperature_2m"])) / 24, d["temperature_2m"], lw=0.9, label=nombre)
    ax[1].axvline(14, color="k", ls=":"); ax[1].text(14.2, ax[1].get_ylim()[1] - 1, "14 d", fontsize=8)
    ax[1].set(xlabel="días desde t0 (jueves 19-11-2026 10:00, UTC−3)", ylabel="temperatura del aire [°C]",
              title="Historias GEFS (NOAA, descarga 16-11-2026, lat −31,1115 lon −71,6147)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "curvas_hidratacion_y_clima.png", dpi=140); plt.close(fig)


def fig_gantt(ids):
    fig, ax = plt.subplots(figsize=(11, 0.5 * len(ids) * 1.1 + 1.5))
    colores = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974"]
    for y, pid in enumerate(ids):
        archivos = sorted(CACHE.glob(f"{pid}__*.json"))
        if not archivos:
            continue
        d = json.loads(archivos[0].read_text(encoding="utf-8"))
        h, sp = d["horarios"], d["plan"]
        for i in range(1, 6):
            a, du = h["activaciones_h"][str(i)], h["duraciones_h"][str(i)]
            ax.barh(y, du, left=a, color=colores[i - 1], edgecolor="k", height=0.55)
            ax.text(a + du / 2, y, f"V{i}", ha="center", va="center", fontsize=7, color="w")
        mold = sp.get("moldaje")
        if mold:
            for i in mold.get("vertidas", []):
                ax.plot([h["activaciones_h"][str(i)] + 72.0], [y - 0.38], marker="D", color="k", ms=4)
        txt = f"D = {h['D_h']:.1f} h"
        if sp.get("manta"):
            ini = sp["manta"].get("inicio_h", 15.0)
            ax.plot([ini], [y + 0.38], marker="v", color="darkorange")
            r = sp["manta"].get("retiro_h")
            if r is None:
                txt += "; manta sin retiro en el horizonte"
            else:
                txt += f"; retiro de manta a las {ini + r:.0f} h"
                if ini + r <= 110:
                    ax.plot([ini + r], [y + 0.38], marker="^", color="crimson")
        ax.text(h["D_h"] + 0.5, y, txt, va="center", fontsize=7)
    ax.set_yticks(range(len(ids))); ax.set_yticklabels(ids, fontsize=8)
    ax.set(xlabel="horas desde t0 (jueves 19-11-2026 10:00)", xlim=(0, 110),
           title="Cronogramas: colocación (barras), instalación de manta ▼, retiro de moldaje ◆ (72 h desde colocación), retiro de manta ▲")
    ax.grid(axis="x", alpha=0.3); ax.invert_yaxis()
    fig.tight_layout(); fig.savefig(FIG / "cronogramas.png", dpi=140); plt.close(fig)


def fig_series(ids, curva, clima):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for pid in ids:
        d = corrida(pid, curva, clima)
        s = d["serie"]
        t = np.asarray(s["h"]) / 24
        ax[0].plot(t, s["Tmax_c"], label=f"{pid}: máx. hormigón")
        ax[0].plot(t, s["S25_min_c"], "--", label=f"{pid}: mín. S25")
        ax[1].plot(t, s["delta_c"], label=pid)
    ax[0].axhline(70, color="r", ls=":", lw=1.5); ax[0].text(0.2, 70.8, "70 °C", color="r")
    ax[1].axhline(20, color="r", ls=":", lw=1.5); ax[1].text(0.2, 20.4, "20 °C", color="r")
    for a in ax:
        a.axvline(14, color="gray", ls=":"); a.axvline(28, color="gray", ls=":"); a.grid(alpha=0.3)
        a.set_xlabel("días desde t0")
    ax[0].set(ylabel="T [°C]", title=f"Temperatura máxima y mínima S25 — {curva}, clima {clima}")
    ax[1].set(ylabel="ΔT = Tmax − T_S25,min [°C]", title="Diferencial núcleo–superficie (S25)")
    ax[0].legend(fontsize=6); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / f"series_{curva}_{clima}.png", dpi=140); plt.close(fig)


def fig_pareto(ruta):
    filas = list(csv.DictReader(Path(ruta).open(encoding="utf-8")))
    interps = sorted({f["interpretacion"] for f in filas})
    fig, axs = plt.subplots(1, len(interps), figsize=(5.2 * len(interps), 4.6), squeeze=False)
    for ax, ip in zip(axs[0], interps):
        for f in filas:
            if f["interpretacion"] != ip:
                continue
            D = float(f["D_h"])
            ax.scatter(D, 100 * float(f["P_f_14d"]), marker="o", color="#4c72b0", s=45)
            k = "P_f_fin" if "P_f_fin" in f and f["P_f_fin"] != "" else "P_f_28d"
            ax.scatter(D, 100 * float(f[k]), marker="s", color="#c44e52", s=45)
            ax.annotate(f["plan_id"].replace("f_", ""), (D, 100 * float(f[k])), fontsize=5, rotation=25)
        ax.axhline(5, color="k", ls="--"); ax.text(0.5, 6, "5 %", fontsize=8)
        ax.set(xlabel="D [h] (fin de colocación − t0)", ylabel="P_f conjunta [%]", title=ip)
        ax.grid(alpha=0.3)
    axs[0][0].scatter([], [], marker="o", color="#4c72b0", label="horizonte 14 d")
    axs[0][0].scatter([], [], marker="s", color="#c44e52", label="horizonte extendido")
    axs[0][0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "D_vs_Pf.png", dpi=140); plt.close(fig)


def fig_cortes(ruta_npz, salida, vmin=10.0, vmax=60.0):
    d = np.load(ruta_npz)
    claves = sorted(k for k in d.files if k.startswith("T_"))
    vid = d["vertida_id"]
    origen, paso = d["origen"], d["paso"]
    fig, axs = plt.subplots(len(claves), 3, figsize=(15, 3.6 * len(claves)), squeeze=False,
                            constrained_layout=True, gridspec_kw={"width_ratios": [1.6, 1.6, 1.0]})
    for fila, k in enumerate(claves):
        T = np.where(d["activa_" + k[2:]] & (vid > 0), d[k], np.nan)
        nx, ny, nz = T.shape
        ix, iy = nx // 2, ny // 2
        iz = int(round((1.65 - origen[2]) / paso[2] - 0.5))
        ext_x = [origen[0], origen[0] + nx * paso[0]]; ext_y = [origen[1], origen[1] + ny * paso[1]]
        ext_z = [origen[2], origen[2] + nz * paso[2]]
        cortes = [(T[:, iy, :].T, ext_x + ext_z, "corte vertical y = 0 (x–z)"),
                  (T[ix, :, :].T, ext_y + ext_z, "corte vertical x = 0 (y–z)"),
                  (T[:, :, iz].T, ext_x + ext_y, "corte horizontal z ≈ 1,65 m (x–y)")]
        for col, (a, ext, tit) in enumerate(cortes):
            im = axs[fila][col].imshow(a, origin="lower", extent=ext, vmin=vmin, vmax=vmax, cmap="inferno",
                                       aspect="equal")
            axs[fila][col].set_title(f"{k[2:]} — {tit}", fontsize=8)
    fig.colorbar(im, ax=axs, label="T [°C] (misma escala 10–60 °C en todos los cortes)", shrink=0.9)
    fig.suptitle(Path(ruta_npz).name.split("__")[0] + " — cortes del campo 3D de temperatura", fontsize=10)
    fig.savefig(salida, dpi=130); plt.close(fig)


def main() -> int:
    modo = sys.argv[1]
    if modo == "base":
        fig_base()
    elif modo == "gantt":
        fig_gantt(sys.argv[2].split(","))
    elif modo == "series":
        fig_series([sys.argv[2]] + sys.argv[5:], sys.argv[3], sys.argv[4])
    elif modo == "pareto":
        fig_pareto(sys.argv[2])
    elif modo == "cortes":
        fig_cortes(sys.argv[2], sys.argv[3])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
