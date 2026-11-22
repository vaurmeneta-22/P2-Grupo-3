"""Compara espesores del Plan 1 con manta mantenida hasta la hora 672."""
from fractions import Fraction
import gzip
import json

import matplotlib.pyplot as plt

from cierre_plan1 import ROOT, SALIDA, SPEC, TABLAS, bateria
from resumir_cierre_plan1 import fila_resumen, guardar_csv

HORIZONTE_H = 672
ESPESORES = (
    (60, "pur_mantenido", None),
    (80, "pur_80mm", None),
    (100, "pur_100mm", "pur_100mm_fino"),
    (120, "pur_120mm", "pur_120mm_mantenido_1008_fino"),
)


def leer_matriz(variante, matriz, dx):
    rutas = {p.name.removesuffix(".json.gz"): p for p in (SALIDA / variante).glob("*.json.gz")}
    if len(matriz) != 24 or rutas.keys() != matriz.keys():
        raise ValueError(f"{variante}: se requieren los 24 escenarios, hay {len(rutas)}")
    filas = {}
    for escenario, ruta in sorted(rutas.items()):
        with gzip.open(ruta, "rt", encoding="utf-8") as archivo:
            datos = json.load(archivo)
        fila = datos["fila"].copy()
        if (fila["escenario"] != escenario or fila["dx_m"] != dx
                or fila["dt_limite_s"] != 600 or fila["peso"] != matriz[escenario][3]
                or fila["horizonte_h"] < HORIZONTE_H or not 0 < fila["dt_sobre_estable_max"] < 1
                or fila["balance_relativo"] >= 1e-12):
            raise ValueError(f"Resultados incompatibles: {ruta}")
        if any(t <= HORIZONTE_H * 3600 for t, nombres, *_ in datos["retiros"]
               if any("manta_PUR" in nombre for nombre in nombres)):
            raise ValueError(f"Retiro de PUR dentro de la ventana: {ruta}")
        serie = [m for m in datos["serie"] if m["tiempo_h"] <= HORIZONTE_H]
        if serie[-1]["tiempo_h"] != HORIZONTE_H:
            raise ValueError(f"Falta la muestra final: {ruta}")
        pico_t = max(serie, key=lambda m: m["temperatura_maxima_c"])
        pico_delta = max(serie, key=lambda m: m["diferencia_maxima_S25_c"])
        t, delta = pico_t["temperatura_maxima_c"], pico_delta["diferencia_maxima_S25_c"]
        if fila["horizonte_h"] == HORIZONTE_H and (t, delta) != (fila["M_T_c"], fila["M_delta_c"]):
            raise ValueError(f"Los extremos no coinciden con la serie: {ruta}")
        fila.update(M_T_c=t, t_M_T_h=pico_t["tiempo_h"], M_delta_c=delta,
                    t_M_delta_h=pico_delta["tiempo_h"], falla_T=t > 70, falla_delta=delta > 20,
                    falla_conjunta=t > 70 or delta > 20)
        filas[escenario] = fila
    return filas


def comprobar_clima(contexto, variante):
    """Verifica que el tramo reutilizado conserva las mismas muestras GEFS."""
    original = ROOT / SPEC["climas_extendidos"]
    archivo = ROOT / SPEC["variantes"][variante]["climas_extendidos"]
    for config in sorted(original.glob("*.jsonc")):
        pronosticos = []
        for ruta in (config, archivo / config.name):
            clima = contexto["ev"].leer_jsonc(ruta)["clima"]
            forecast = ROOT / clima["snapshot"]["carpeta"] / "forecast.json"
            pronosticos.append(json.loads(forecast.read_text()))
        base, extendido = pronosticos
        for clave in ("time_utc", "temperature_2m"):
            if len(base[clave]) != HORIZONTE_H + 1 or base[clave] != extendido[clave][:HORIZONTE_H + 1]:
                raise ValueError(f"El clima a 28 días difiere en {config.name}: {clave}")


def comprobar_recorte():
    escenario = "M4_SS__GEFS_20261116_central.json.gz"
    series = []
    for variante in ("pur_120mm_fino", "pur_120mm_mantenido_1008_fino"):
        with gzip.open(SALIDA / variante / escenario, "rt", encoding="utf-8") as archivo:
            serie = json.load(archivo)["serie"]
        series.append([m for m in serie if m["tiempo_h"] <= HORIZONTE_H])
    directa, recortada = series
    if directa != recortada[:len(directa)] or any(m != directa[-1] for m in recortada[len(directa):]):
        raise ValueError("La serie de control de 120 mm difiere al recortar la extensión")


def graficar(resumen):
    fig, eje = plt.subplots(figsize=(8, 4.5), layout="constrained")
    for dx, etiqueta, marcador in ((0.25, "Malla 0,25 m", "o"), (0.125, "Malla 0,125 m", "s")):
        filas = [fila for fila in resumen if fila["dx_m"] == dx]
        eje.plot([fila["espesor_mm"] for fila in filas], [fila["max_delta_c"] for fila in filas],
                 linestyle="none", marker=marcador, label=etiqueta)
    eje.axhline(20, color="0.35", linestyle="--", label="Umbral por escenario")
    eje.set(xlabel="Espesor de PUR mantenido [mm]", ylabel="Mayor diferencial S25 [°C]",
            title="Plan 1 a 28 días")
    eje.set_xticks([mm for mm, _, _ in ESPESORES])
    eje.legend()
    fig.savefig(TABLAS / "pur_28d_diferencial.png", dpi=160)
    plt.close(fig)


def main():
    contexto = bateria._contexto()
    matriz = contexto["matriz"]
    comprobar_clima(contexto, "pur_120mm_mantenido_1008_fino")
    comprobar_recorte()
    resumen, comparacion = [], []
    for mm, variante_nominal, variante_fina in ESPESORES:
        por_malla = {}
        for variante, dx in ((variante_nominal, 0.25), (variante_fina, 0.125)):
            if variante is None:
                continue
            filas = leer_matriz(variante, matriz, dx)
            por_malla[dx] = filas
            resumen.append(dict(espesor_mm=mm, dx_m=dx, corte_h=HORIZONTE_H,
                                horizonte_corrida_h=next(iter(filas.values()))["horizonte_h"],
                                **fila_resumen(variante, list(filas.values()), matriz)))
        if variante_fina is None:
            continue
        for escenario in sorted(matriz):
            peso = Fraction(matriz[escenario][3]).limit_denominator(248)
            fila = dict(espesor_mm=mm, escenario=escenario, peso=str(peso))
            for dx, sufijo in ((0.25, "025"), (0.125, "0125")):
                datos = por_malla[dx][escenario]
                for clave in ("M_T_c", "t_M_T_h", "M_delta_c", "t_M_delta_h",
                              "falla_T", "falla_delta", "falla_conjunta"):
                    fila[f"{clave}_{sufijo}"] = datos[clave]
            comparacion.append(fila)
    guardar_csv("pur_28d_resumen.csv", resumen)
    graficar(resumen)
    for mm, _, variante_fina in ESPESORES:
        if variante_fina:
            guardar_csv(f"pur_28d_{mm}mm_por_escenario.csv",
                        [fila for fila in comparacion if fila["espesor_mm"] == mm])


if __name__ == "__main__":
    main()
