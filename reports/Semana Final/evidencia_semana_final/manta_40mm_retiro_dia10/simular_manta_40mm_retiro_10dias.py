from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(r"C:\Users\56990\Desktop\Universidad\8vo Semestre\Metodos Computacionales\P2-Grupo-3")
OUT = Path(r"C:\Users\56990\Documents\Codex\2026-10-15\h\outputs\recalculo_40mm_retiro_10dias")
BASE = Path(r"C:\Users\56990\Documents\Codex\2026-10-15\h\outputs\recalculo_40mm_retiro_10dias\configuracion_usada.json")
CFG = BASE
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import demo_e10_entorno as demo
import demo_pipeline_ifc_3d as pipe
import evaluar_dossier_final as ev
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.decision import SeguimientoMetricasDecision
from mcplan.weather.loader import construir_aire

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def shifted_neighbor(mask, axis, direction):
    shifted = np.roll(mask, direction, axis=axis)
    sl = [slice(None)] * mask.ndim
    sl[axis] = 0 if direction == 1 else -1
    shifted[tuple(sl)] = False
    return shifted

class SerieGlobal:
    def __init__(self, grilla, profundidad):
        self.detector = SeguimientoMetricasDecision(grilla, profundidad)
        self.grilla = grilla
        self.rows = []

    def registrar(self, t, T, solver):
        self.detector.registrar(t, T, solver)
        active = solver.mascara_activa.reshape(-1)
        field = np.asarray(T, dtype=float).reshape(-1)
        ids = np.flatnonzero(active)
        if not ids.size:
            return
        i_hot = int(ids[np.argmax(field[ids])])
        hot_xyz = np.asarray(self.grilla.origen_m) + (
            np.asarray(np.unravel_index(i_hot, self.grilla.forma)) + 0.5
        ) * np.asarray(self.grilla.paso_m)

        s25_values, s25_xyz, s25_faces = [], [], []
        for face, indices0, condition, rc0, rx0, d, xyz0 in self.detector._partes:
            keep = active[indices0]
            if not np.any(keep):
                continue
            indices = indices0[keep]
            rc = np.asarray(rc0)[keep]
            rx = np.asarray(rx0)[keep]
            tc = field[indices]
            tf = self.detector._cara_c(tc, condition, rc, rx, float(t))
            s25 = tf + (self.detector.profundidad_m / (d / 2.0)) * (tc - tf)
            s25_values.append(s25)
            s25_xyz.append(xyz0[keep])
            s25_faces.extend([face] * int(np.count_nonzero(keep)))
        if not s25_values:
            return
        values = np.concatenate(s25_values)
        xyz = np.concatenate(s25_xyz)
        j = int(np.argmin(values))
        row = {
            "t_h": float(t) / 3600.0,
            "Tmax_c": float(field[ids].max()),
            "Tmin_c": float(field[ids].min()),
            "Tmedia_c": float(field[ids].mean()),
            "S25_min_c": float(values[j]),
            "DeltaT_global_c": float(field[i_hot] - values[j]),
            "xyz_caliente_m": [float(v) for v in hot_xyz],
            "T_caliente_c": float(field[i_hot]),
            "xyz_S25_frio_m": [float(v) for v in xyz[j]],
            "cara_S25_frio": s25_faces[j],
        }
        # El solver registra un estado antes y después de un evento en el mismo instante.
        if self.rows and abs(self.rows[-1]["t_h"] - row["t_h"]) < 1e-12:
            self.rows[-1] = row
        else:
            self.rows.append(row)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.loads(BASE.read_text(encoding="utf-8"))
    try:
        config, dossier = ev.cargar_configuracion_dossier(CFG)
        clima_cfg = dossier["climas"][0]
        overlay = ev.leer_jsonc((ROOT / clima_cfg["config_jsonc"]).resolve())
        clima = dict(overlay.get("clima") or {})
        t0 = dossier["inicio_obra_utc"]
        clima["inicio_obra_utc"] = t0
        duration_s = float(dossier["tiempo_final_h"]) * 3600.0
        air_check = construir_aire(clima, tiempo_final_s=duration_s, raiz_repo=ROOT)
        if air_check.es_sintetica:
            raise RuntimeError("El caso debe usar el snapshot de clima GEFS central.")
        config_run = dict(config)
        for key in ("clima", "sitio", "entorno", "condicion_superficie"):
            if key in overlay:
                config_run[key] = overlay[key]
        config_run.setdefault("clima", {})["inicio_obra_utc"] = t0

        originals = extraer_vertidas((ROOT / config_run["ifc"]).resolve())
        plan = ev.cargar_plan(dossier["planes"][0], originals)
        pours = plan.aplicar(originals)
        grid = voxelizar_vertidas(pours, float(dossier["paso_voxel_m"]))
        observer = SerieGlobal(grid, float(dossier.get("profundidad_s25_m", 0.025)))
        result, borders, air = demo.correr(
            config_run, grid, pours, observadores=(observer,),
            modelo_hidratacion=ev.leer_curvas(dossier)[0][1])

        rows = observer.rows
        if not rows:
            raise RuntimeError("El solver no registró estados térmicos.")
        delta_row = max(rows, key=lambda r: r["DeltaT_global_c"])
        tmax_row = max(rows, key=lambda r: r["Tmax_c"])
        avg_row = max(rows, key=lambda r: r["Tmedia_c"])
        metricas_oficiales = observer.detector.resumen()
        if not np.isclose(delta_row["DeltaT_global_c"], metricas_oficiales["M_delta_c"], atol=1e-9):
            raise AssertionError("La serie de ΔT no coincide con SeguimientoMetricasDecision.")
        if not np.isclose(tmax_row["Tmax_c"], metricas_oficiales["M_T_c"], atol=1e-9):
            raise AssertionError("La serie de Tmax no coincide con SeguimientoMetricasDecision.")
        if not np.isclose(rows[-1]["t_h"], float(dossier["tiempo_final_h"]), atol=1e-9):
            raise AssertionError("La serie no alcanza el final declarado de la simulación.")
        with (OUT / "series_temperatura.csv").open("w", newline="", encoding="utf-8") as f:
            fields = ["t_h", "Tmax_c", "Tmin_c", "Tmedia_c", "S25_min_c", "DeltaT_global_c",
                      "cara_S25_frio", "x_caliente_m", "y_caliente_m", "z_caliente_m",
                      "x_S25_frio_m", "y_S25_frio_m", "z_S25_frio_m"]
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for r in rows:
                hot, cold = r["xyz_caliente_m"], r["xyz_S25_frio_m"]
                writer.writerow({
                    "t_h": r["t_h"], "Tmax_c": r["Tmax_c"], "Tmin_c": r["Tmin_c"],
                    "Tmedia_c": r["Tmedia_c"],
                    "S25_min_c": r["S25_min_c"], "DeltaT_global_c": r["DeltaT_global_c"],
                    "cara_S25_frio": r["cara_S25_frio"],
                    "x_caliente_m": hot[0], "y_caliente_m": hot[1], "z_caliente_m": hot[2],
                    "x_S25_frio_m": cold[0], "y_S25_frio_m": cold[1], "z_S25_frio_m": cold[2]})

        horas = np.array([r["t_h"] for r in rows])
        tmax = np.array([r["Tmax_c"] for r in rows])
        tmin = np.array([r["Tmin_c"] for r in rows])
        tmean = np.array([r["Tmedia_c"] for r in rows])
        s25 = np.array([r["S25_min_c"] for r in rows])
        delta = np.array([r["DeltaT_global_c"] for r in rows])
        plt.rcParams.update({"font.size": 10, "axes.titleweight": "semibold"})
        fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 1.0]})
        axes[0].plot(horas, tmax, color="#d95f02", lw=2.2, label="T máxima")
        axes[0].plot(horas, tmean, color="#1b9e77", lw=2.0, label="T media del hormigón activo")
        axes[0].plot(horas, tmin, color="#7570b3", lw=1.8, label="T mínima del hormigón activo")
        axes[0].plot(horas, s25, color="#377eb8", lw=1.6, ls="--", label="S25 mínima global")
        axes[0].scatter([tmax_row["t_h"]], [tmax_row["Tmax_c"]], color="#d95f02", zorder=4)
        axes[0].set_ylabel("Temperatura [°C]")
        axes[0].set_title("Evolución de temperaturas · manta PUR 40 mm · retiro a 10 días")
        axes[0].grid(True, alpha=0.25)
        axes[0].legend(ncol=2, frameon=True)
        axes[1].plot(horas, delta, color="#6a3d9a", lw=2.1, label="ΔT global = Tmax − S25 mínima")
        axes[1].axhline(20.0, color="#b2182b", ls="--", lw=1.4, label="Referencia del proyecto: 20 °C")
        axes[1].scatter([delta_row["t_h"]], [delta_row["DeltaT_global_c"]], color="#6a3d9a", zorder=4)
        axes[1].annotate(f'Máx. {delta_row["DeltaT_global_c"]:.2f} °C · {delta_row["t_h"]:.1f} h',
                         (delta_row["t_h"], delta_row["DeltaT_global_c"]),
                         xytext=(-125, 12), textcoords="offset points",
                         arrowprops={"arrowstyle":"->", "color":"#6a3d9a"})
        axes[1].set_xlabel("Tiempo desde el inicio de obra [h]")
        axes[1].set_ylabel("Diferencia [°C]")
        axes[1].grid(True, alpha=0.25)
        axes[1].legend(frameon=True)
        fig.tight_layout()
        fig.savefig(OUT / "evolucion_temperatura_40mm_retiro_10dias.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

        active = np.asarray(grid.vertida_id) > 0
        all_neighbors = np.ones_like(active)
        for axis in range(3):
            all_neighbors &= shifted_neighbor(active, axis, 1)
            all_neighbors &= shifted_neighbor(active, axis, -1)
        shell = active & ~all_neighbors
        ijk = np.argwhere(shell)
        xyz_shell = np.asarray(grid.origen_m) + (ijk + 0.5) * np.asarray(grid.paso_m)
        hot = np.asarray(delta_row["xyz_caliente_m"])
        cold = np.asarray(delta_row["xyz_S25_frio_m"])
        fig = plt.figure(figsize=(14, 7.5))
        ax3 = fig.add_subplot(1, 2, 1, projection="3d")
        ax2 = fig.add_subplot(1, 2, 2)
        ax3.scatter(xyz_shell[:,0], xyz_shell[:,1], xyz_shell[:,2], s=3, c="#9aa7b1",
                    alpha=0.22, depthshade=False, label="Superficie voxelizada")
        ax3.scatter(*hot, s=150, marker="*", c="#e4572e", edgecolors="black",
                    linewidths=0.6, label="Más caliente en el instante del ΔT máximo")
        ax3.scatter(*cold, s=120, marker="*", c="#2878b5", edgecolors="black",
                    linewidths=0.6, label="S25 más frío")
        ax3.plot([hot[0], cold[0]], [hot[1], cold[1]], [hot[2], cold[2]],
                 color="#745296", alpha=0.7, lw=1.4, ls=":")
        ax3.set_xlabel("X [m]"); ax3.set_ylabel("Y [m]"); ax3.set_zlabel("Z [m]")
        ax3.set_title(f"Ubicación de los extremos · t = {delta_row['t_h']:.1f} h")
        ax3.view_init(elev=25, azim=-55)
        ax3.legend(loc="upper left", fontsize=8)
        ax2.scatter(xyz_shell[:,0], xyz_shell[:,1], s=3, c="#b5bec5", alpha=0.25,
                    label="Proyección de la superficie")
        ax2.scatter(hot[0], hot[1], s=150, marker="*", c="#e4572e", edgecolors="black",
                    linewidths=0.6, label=f"Caliente ({hot[0]:.2f}, {hot[1]:.2f}, {hot[2]:.2f}) m")
        ax2.scatter(cold[0], cold[1], s=140, marker="*", c="#2878b5", edgecolors="black",
                    linewidths=0.6, label=f"S25 frío ({cold[0]:.2f}, {cold[1]:.2f}, {cold[2]:.2f}) m")
        ax2.set_aspect("equal", adjustable="box")
        ax2.set_xlabel("X [m]"); ax2.set_ylabel("Y [m]")
        ax2.set_title("Vista en planta (proyección XY)")
        ax2.grid(True, alpha=0.25)
        ax2.legend(loc="best", fontsize=8)
        fig.suptitle(f"ΔT global máximo = {delta_row['DeltaT_global_c']:.2f} °C · "
                     f"Tmax = {delta_row['T_caliente_c']:.2f} °C · "
                     f"S25 = {delta_row['S25_min_c']:.2f} °C",
                     fontweight="semibold")
        fig.tight_layout()
        fig.savefig(OUT / "ubicacion_extremos_40mm_retiro_10dias.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

        resumen = {
            "manta_espesor_m": 0.04,
            "retiro_manta_h": 240.0,
            "inicio_obra_local": dossier["inicio_obra_local"],
            "horizonte_h": dossier["tiempo_final_h"],
            "paso_voxel_m": dossier["paso_voxel_m"],
            "paso_tiempo_max_s": dossier["paso_tiempo_max_s"],
            "escenario": f'{dossier["curvas"][0]["id"]}__{clima_cfg["id"]}',
            "numero_muestras_solver": len(rows),
            "validacion_series_vs_observador_oficial": "coinciden; maximos y tiempo final comprobados",
            "DeltaT_global_max_c": delta_row["DeltaT_global_c"],
            "t_DeltaT_global_max_h": delta_row["t_h"],
            "Tmax_c": tmax_row["Tmax_c"],
            "t_Tmax_h": tmax_row["t_h"],
            "Tmin_global_c": float(tmin.min()),
            "t_Tmin_global_h": float(horas[int(np.argmin(tmin))]),
            "Tmedia_max_c": avg_row["Tmedia_c"],
            "t_Tmedia_max_h": avg_row["t_h"],
            "pareja_en_max_DeltaT": {
                "caliente_c": delta_row["T_caliente_c"],
                "caliente_xyz_m": delta_row["xyz_caliente_m"],
                "S25_min_c": delta_row["S25_min_c"],
                "S25_xyz_m": delta_row["xyz_S25_frio_m"],
                "cara_S25": delta_row["cara_S25_frio"],
            },
            "aire": air.descripcion,
            "salidas": ["series_temperatura.csv", "evolucion_temperatura_40mm_retiro_10dias.png",
                        "ubicacion_extremos_40mm_retiro_10dias.png"],
        }
        (OUT / "resumen.json").write_text(json.dumps(resumen, ensure_ascii=False, indent=2)+"\n",
                                          encoding="utf-8")
        print(json.dumps(resumen, ensure_ascii=False, indent=2))
        print(f"Resultados: {OUT}")
    finally:
        pass

if __name__ == "__main__":
    main()



