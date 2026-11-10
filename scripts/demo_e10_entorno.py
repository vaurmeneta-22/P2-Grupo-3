"""Demo corta E10: aire variable, suelo, moldaje/aislación y su retiro.

Ejecutar con ▶️ (o ``python scripts/demo_e10_entorno.py``). Tarda unos segundos.

Qué muestra (figura ``results/e10_demo/demo_e10.png`` y texto en consola):
  1. Temperatura ambiente (SERIE SINTÉTICA), temperatura interior y superficial (cara
     estimada del hormigón) de la vertida 2.
  2. Flujo por la cara lateral x+ de las vertidas 1 y 2 en el tiempo: salta cuando se retira
     la protección (moldaje a las 72 h desde la colocación; aislación de la V1 a las 120 h
     desde el inicio de obra).
  3. Tmax con la serie de aire base frente a la misma serie +10 °C: cambiar el aire
     cambia el resultado. Además comprueba que ``entorno.activo = false`` sigue
     reproduciendo la corrida E09 (mismo resultado que condicion_superficie).

La serie de aire es sintética y las coordenadas del Grupo 3 siguen pendientes de sitios.yaml:
la demo no consulta ninguna API.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
sys.path.insert(0, str(AQUI))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import demo_pipeline_ifc_3d as pipe
from mcplan.geometry.ifc import extraer_vertidas
from mcplan.geometry.voxelizacion import voxelizar_vertidas
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.seguimiento import SeguimientoPuntos, proponer_puntos, resolver_puntos
from mcplan.thermal.vertidas import simular_vertidas

SALIDA = RAIZ / "results" / "e10_demo"
CONFIG_BASE = RAIZ / "configs" / "semana3_fundacion.jsonc"
CONFIG_DEMO = RAIZ / "configs" / "e10_demo.jsonc"


class FlujoLateral:
    """Observa el flujo [W/m²] por la cara x+ de cada vertida (+ = entra al hormigón)."""

    def __init__(self, grilla, ids, cara="x+"):
        self.grilla, self.ids, self.cara = grilla, ids, cara
        self.tiempos_s: list[float] = []
        self.flujo = {i: [] for i in ids}

    def registrar(self, tiempo_s, temperatura, solver):
        if self.tiempos_s and tiempo_s == self.tiempos_s[-1]:
            self.tiempos_s.pop()
            for serie in self.flujo.values():
                serie.pop()
        self.tiempos_s.append(float(tiempo_s))
        flujo_w = solver.flujos_frontera_w(temperatura, tiempo_s)[self.cara]
        expuesta = solver.superficie_expuesta[self.cara]
        area = solver.area_cara_m2["xyz".index(self.cara[0])]
        for i in self.ids:
            celdas = expuesta & (self.grilla.vertida_id == i)
            n = int(celdas.sum())
            self.flujo[i].append(float(flujo_w[celdas].sum() / (n * area)) if n else np.nan)


def cargar_config() -> dict:
    config = pipe.cargar_configuracion(CONFIG_BASE)
    superposicion = json.loads(pipe.quitar_comentarios_jsonc(CONFIG_DEMO.read_text(encoding="utf-8")))
    config.update(superposicion)                              # superpone los bloques E10
    return config


def correr(config: dict, grilla, vertidas, *, observadores=()):
    materiales = config["materiales"]
    k, rho, cp = pipe.campos_materiales(grilla, vertidas, materiales)
    hcfg = config["hidratacion"]
    hidratacion = ModeloPulsoCosenoidal(float(hcfg["qmax_w_m3"]),
                                        float(hcfg["duracion_pulso_h"]) * 3600.0)
    final_s = float(config["tiempo_final_h"]) * 3600.0
    bordes, aire = pipe.crear_bordes_entorno(config, grilla, vertidas, final_s)
    resultado = simular_vertidas(
        grilla=grilla, vertidas=vertidas, conductividad_w_m_k=k, densidad_kg_m3=rho,
        calor_especifico_j_kg_k=cp, hidratacion=hidratacion, tiempo_final_s=final_s,
        fraccion_paso_estable=float(config.get("fraccion_paso_estable", 0.8)),
        paso_tiempo_max_s=float(config.get("paso_tiempo_max_s", 600.0)),
        tiempos_muestreo_s=np.arange(0.0, final_s + 1e-9, 3600.0).tolist(),
        fronteras=pipe.crear_fronteras(config), seguimiento=list(observadores) or None,
        bordes=bordes)
    return resultado, bordes, aire


def main() -> int:
    SALIDA.mkdir(parents=True, exist_ok=True)
    config = cargar_config()
    vertidas = pipe.aplicar_overrides(extraer_vertidas((RAIZ / config["ifc"]).resolve()), config)
    grilla = voxelizar_vertidas(vertidas, float(config["paso_voxel_m"]))
    print(f"Sitio: {config['sitio']['nombre']}, {config['sitio']['region']} "
          f"(lat/lon: {config['sitio']['latitud']}/{config['sitio']['longitud']} -> "
          f"{config['sitio']['fuente_coordenadas']})")

    # --- Corrida principal con observadores ---
    puntos = resolver_puntos(grilla, proponer_puntos(grilla, "x+"))
    seguimiento = SeguimientoPuntos(grilla, puntos)
    flujo = FlujoLateral(grilla, [1, 2, 5])
    resultado, bordes, aire = correr(config, grilla, vertidas, observadores=(seguimiento, flujo))
    horas_flujo = np.asarray(flujo.tiempos_s) / 3600.0
    print(f"Aire: {aire.descripcion['aviso']}")
    print("Retiros (t, protecciones, paso estable antes -> despues):")
    for t_s, nombres, antes, despues in resultado.retiros:
        print(f"  t={t_s / 3600:7.2f} h  {', '.join(nombres)}  {antes:7.1f} s -> {despues:7.1f} s")

    # --- Comparaciones ---
    mas_calido = copy.deepcopy(config)
    mas_calido["clima"]["sintetica"]["media_c"] += 10.0
    res_calido, _, _ = correr(mas_calido, grilla, vertidas)
    sin_entorno = copy.deepcopy(config)
    sin_entorno["entorno"] = {"activo": False}
    res_e09, _, _ = correr(sin_entorno, grilla, vertidas)
    dif_aire = float(np.max(np.abs(res_calido.temperaturas_maximas_c - resultado.temperaturas_maximas_c)))
    print("\nMismo aire, cara x+ protegida (V1) frente a expuesta (V5) [W/m2, + entra]:")
    for hora in (30.0, 60.0, 100.0, 140.0):
        n = int(np.argmin(np.abs(horas_flujo - hora)))
        print(f"  t={hora:5.0f} h  V1={flujo.flujo[1][n]:7.1f}  V5={flujo.flujo[5][n]:7.1f}")
    print(f"\nTmax pico con entorno: {resultado.temperaturas_maximas_c.max():.3f} °C | "
          f"aire +10 °C: {res_calido.temperaturas_maximas_c.max():.3f} °C | "
          f"sin entorno (E09, aire constante {config['condicion_superficie']['temperatura_aire_c']:g} °C): "
          f"{res_e09.temperaturas_maximas_c.max():.3f} °C")
    print(f"Cambiar el aire +10 °C modifica Tmax hasta {dif_aire:.3f} °C (debe ser > 0)")
    assert dif_aire > 0.5, "cambiar la serie de aire debía modificar el resultado"

    # --- Figura ---
    horas_obs = np.asarray(seguimiento.tiempos_s) / 3600.0
    punto = next(p for p in puntos if p.vertida_id == 2)
    serie = seguimiento.serie(2)
    aire_h = [aire.funcion(t * 3600.0) for t in horas_obs]
    fig, ejes = plt.subplots(3, 1, figsize=(10, 11), sharex=False)
    ejes[0].plot(horas_obs, aire_h, color="#888888", label="Aire (SERIE SINTÉTICA)")
    ejes[0].plot(horas_obs, serie["interior_c"], color="#d62728", label="Interior V2")
    ejes[0].plot(horas_obs, serie["cara_estimada_c"], color="#2ca02c", linestyle="--",
                 label=f"Superficie V2 (cara {punto.cara} del hormigón, estimada)")
    ejes[0].set_ylabel("Temperatura [°C]")
    ejes[0].set_title("Aire, interior y superficie (V2)")
    horas_flujo = np.asarray(flujo.tiempos_s) / 3600.0
    ejes[1].plot(horas_flujo, flujo.flujo[1], color="#1f77b4", label="V1 (moldaje + aislación)")
    ejes[1].plot(horas_flujo, flujo.flujo[2], color="#ff7f0e", label="V2 (moldaje)")
    ejes[1].plot(horas_flujo, flujo.flujo[5], color="#2ca02c", label="V5 (expuesta, sin protección)")
    ejes[1].axhline(0.0, color="#555555", linewidth=0.8)
    ejes[1].set_ylabel("Flujo por cara x+ [W/m²]\n(+ entra, − sale)")
    ejes[1].set_title("Flujo por la cara lateral: salta al retirar la protección")
    ejes[2].plot(resultado.tiempos_s / 3600.0, resultado.temperaturas_maximas_c, label="Entorno E10 (aire sintético)")
    ejes[2].plot(res_calido.tiempos_s / 3600.0, res_calido.temperaturas_maximas_c, linestyle="--",
                 label="Mismo entorno, aire +10 °C")
    ejes[2].plot(res_e09.tiempos_s / 3600.0, res_e09.temperaturas_maximas_c, linestyle=":",
                 label="Sin entorno (E09: Robin con aire constante)")
    ejes[2].set_ylabel("Tmax [°C]")
    ejes[2].set_title("Cambiar la serie de aire cambia el resultado")
    for eje in ejes:
        for t_s, nombres, _, _ in resultado.retiros:
            eje.axvline(t_s / 3600.0, color="#999999", linestyle=":", linewidth=1.0)
        eje.set_xlabel("Tiempo desde inicio de obra [h]")
        eje.grid(True, alpha=0.3)
        eje.legend(fontsize=8)
    fig.tight_layout()
    ruta = SALIDA / "demo_e10.png"
    fig.savefig(ruta, dpi=140)
    print(f"\nFigura: {ruta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
