"""Diferencias térmicas máximas entre celdas y contra la superficie exterior.

Observa la MISMA simulación (tras cada paso de integración y cada activación) y,
para el dominio activo completo y para cada vertida por separado, calcula:

1. ``ΔT_global = Tmax − Tmin`` entre CENTROS de voxel activos, con las celdas que
   producen ambos extremos. También registra la temperatura media de las celdas
   activas (voxeles de igual volumen: media aritmética).
2. Contraste con la superficie exterior: se reconstruye la temperatura de cada
   cara exterior expuesta (con las resistencias de media celda y las condiciones
   de borde del solver; ver ``seguimiento.temperatura_cara_estimada``), se toman
   los extremos de caras ``Tcara_min`` y ``Tcara_max`` y el mayor contraste
   absoluto celda-cara es

       max( |Tcelda_max − Tcara_min| , |Tcelda_min − Tcara_max| )

   (el máximo de |a − b| sobre dos conjuntos se alcanza en los extremos, así que no
   se forman pares). Se conserva el signo ``Tcelda − Tcara``, las ubicaciones de
   ambos puntos y la orientación de la cara.

Una vertida se compara solo consigo misma (sus celdas y sus caras expuestas); si no
está activa, o no tiene caras expuestas, el resultado es n/a (nunca 0 °C). Las
celdas inactivas se excluyen de todos los extremos. Una cara que pasa a ser
interfaz con otra vertida activa deja de ser exterior y sale de los conjuntos de
caras desde ese instante.

ADVERTENCIAS DE INTERPRETACIÓN
* ΔT_global compara temperaturas de centros de celda, no de superficies.
* El contraste superficial usa una temperatura de cara RECONSTRUIDA, no medida.
* Los dos puntos de un extremo pueden estar lejos, no ser vecinos ni estar alineados.
* Son diferencias en °C (o K) entre puntos, no gradientes en K/m, y no son un
  criterio de fisuración ni un análisis mecánico.
* No guarda campos completos: por instante solo conserva escalares e índices.
  Las posiciones de los extremos pueden cambiar con el tiempo (a diferencia de los
  puntos fijos de ``seguimiento``).
* "Todos los pasos" (cada paso de integración y cada activación) y "fotogramas"
  (solo los instantes de un subconjunto, p. ej. los cuadros del GIF) pueden dar
  máximos distintos: el máximo entre fotogramas nunca supera al de todos los pasos.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from .seguimiento import CARAS, NA, _centro_voxel

_CAMPOS = ("tmin", "cel_min", "tmax", "cel_max", "n_celdas",
           "fmin", "fmin_cel", "fmin_cara", "fmax", "fmax_cel", "fmax_cara", "n_caras",
           "con_signo", "con_cel", "con_tcel", "con_fcel", "con_fcara", "con_tcara", "tmedia")
_VACIA = (math.nan, -1, math.nan, -1, 0, math.nan, -1, -1, math.nan, -1, -1, 0,
          math.nan, -1, math.nan, -1, -1, math.nan, math.nan)


def _fmt(valor, decimales=6):
    return NA if valor is None or not math.isfinite(valor) else f"{valor:.{decimales}f}"


class SeguimientoExtremos:
    """Registra extremos de celdas y contraste con caras exteriores por alcance."""

    def __init__(self, grilla):
        self.grilla = grilla
        self.ids = sorted(int(i) for i in np.unique(grilla.vertida_id) if i > 0)
        self.alcances = ["global"] + [f"vertida_{i}" for i in self.ids]
        self.tiempos_s: list[float] = []
        self._filas = {a: [] for a in self.alcances}
        self._ids_planos = grilla.vertida_id.reshape(-1)
        self._referencia = None            # identifica la máscara/exposición ya preparada

    # -- preparación cacheada (solo cambia con cada activación) ---------------
    def _preparar(self, solver):
        self._referencia = solver.superficie_expuesta["x+"]
        activas = solver.mascara_activa.reshape(-1)
        self._celdas = {"global": np.flatnonzero(activas)}
        for i in self.ids:
            self._celdas[f"vertida_{i}"] = np.flatnonzero(activas & (self._ids_planos == i))
        celda, codigo, resistencia, extra, segmentos, inicio = [], [], [], [], [], 0
        for numero, cara in enumerate(CARAS):
            eje = "xyz".index(cara[0])
            r_extra = solver.resistencias_borde.get(cara)
            for seleccion, condicion in solver.tramos_borde[cara]:
                indices = np.flatnonzero(seleccion)
                resistencia.append(solver.dx_m[eje] / (2.0 * solver.k.reshape(-1)[indices]))
                extra.append(np.zeros(len(indices)) if r_extra is None
                             else np.asarray(r_extra).reshape(-1)[indices])
                celda.append(indices)
                codigo.append(np.full(len(indices), numero, dtype=np.int8))
                segmentos.append((slice(inicio, inicio + len(indices)), condicion))
                inicio += len(indices)
        self._cel_cara = np.concatenate(celda)
        self._cod_cara = np.concatenate(codigo)
        self._resistencia = np.concatenate(resistencia)
        self._extra = np.concatenate(extra)
        self._segmentos = segmentos
        vid_cara = self._ids_planos[self._cel_cara]
        self._entradas = {"global": np.arange(len(self._cel_cara))}
        for i in self.ids:
            self._entradas[f"vertida_{i}"] = np.flatnonzero(vid_cara == i)

    def _temperaturas_cara(self, temperatura_plana, tiempo_s):
        t_celda = temperatura_plana[self._cel_cara]
        salida = np.empty_like(t_celda)
        for tramo, cond in self._segmentos:
            tc, rc, rx = t_celda[tramo], self._resistencia[tramo], self._extra[tramo]
            if cond.tipo == "dirichlet":
                salida[tramo] = cond.temperatura(tiempo_s)
            elif cond.tipo == "neumann":
                salida[tramo] = tc + cond.flujo_entrante(tiempo_s) * rc
            elif cond.coeficiente_h_w_m2_k <= 0.0:           # Robin con h = 0: aislamiento
                salida[tramo] = tc
            else:
                salida[tramo] = tc + (cond.temperatura(tiempo_s) - tc) * rc / (
                    rc + rx + 1.0 / cond.coeficiente_h_w_m2_k)
        return salida

    # -- registro ---------------------------------------------------------------
    def registrar(self, tiempo_s: float, temperatura_c: np.ndarray, solver) -> None:
        """Toma una muestra; si repite el último instante la reemplaza (activaciones)."""
        if self._referencia is not solver.superficie_expuesta["x+"]:
            self._preparar(solver)
        if self.tiempos_s and tiempo_s == self.tiempos_s[-1]:
            self.tiempos_s.pop()
            for filas in self._filas.values():
                filas.pop()
        self.tiempos_s.append(float(tiempo_s))
        plana = temperatura_c.reshape(-1)
        t_cara = self._temperaturas_cara(plana, tiempo_s)
        for alcance in self.alcances:
            celdas, entradas = self._celdas[alcance], self._entradas[alcance]
            if len(celdas) == 0:
                self._filas[alcance].append(_VACIA)
                continue
            t_celdas = plana[celdas]
            i_min, i_max = int(np.argmin(t_celdas)), int(np.argmax(t_celdas))
            tmin, cel_min = float(t_celdas[i_min]), int(celdas[i_min])
            tmax, cel_max = float(t_celdas[i_max]), int(celdas[i_max])
            valores = _VACIA[5:-1]
            if len(entradas):
                t_caras = t_cara[entradas]
                j_min, j_max = int(np.argmin(t_caras)), int(np.argmax(t_caras))
                e_min, e_max = entradas[j_min], entradas[j_max]
                a = tmax - float(t_caras[j_min])      # Tcelda_max − Tcara_min
                b = tmin - float(t_caras[j_max])      # Tcelda_min − Tcara_max
                if abs(a) >= abs(b):
                    contraste = (a, cel_max, tmax, e_min, float(t_caras[j_min]))
                else:
                    contraste = (b, cel_min, tmin, e_max, float(t_caras[j_max]))
                valores = (float(t_caras[j_min]), int(self._cel_cara[e_min]), int(self._cod_cara[e_min]),
                           float(t_caras[j_max]), int(self._cel_cara[e_max]), int(self._cod_cara[e_max]),
                           len(entradas), contraste[0], contraste[1], contraste[2],
                           int(self._cel_cara[contraste[3]]), int(self._cod_cara[contraste[3]]),
                           contraste[4])
            self._filas[alcance].append((tmin, cel_min, tmax, cel_max, len(celdas), *valores,
                                         float(np.mean(t_celdas))))

    # -- consulta -----------------------------------------------------------------
    def serie(self, alcance: str) -> dict:
        """Series numpy por campo (NaN = n/a) más ``dt_global_c`` y ``contraste_abs_c``."""
        datos = np.array(self._filas[alcance], dtype=float).reshape(-1, len(_CAMPOS))
        serie = {campo: datos[:, n] for n, campo in enumerate(_CAMPOS)}
        serie["tiempos_s"] = np.asarray(self.tiempos_s, dtype=float)
        serie["dt_global_c"] = serie["tmax"] - serie["tmin"]
        serie["contraste_abs_c"] = np.abs(serie["con_signo"])
        return serie

    def _coord(self, plano):
        return _centro_voxel(self.grilla, np.unravel_index(int(plano), self.grilla.forma))

    def _coord_cara(self, plano, codigo):
        centro = list(self._coord(plano))
        cara = CARAS[int(codigo)]
        eje = "xyz".index(cara[0])
        centro[eje] += (-0.5 if cara[1] == "-" else 0.5) * self.grilla.paso_m[eje]
        return tuple(centro), cara

    def resumen(self, tiempos_fotogramas_s=None) -> list[dict]:
        """Máximos temporales por alcance y métrica, en todos los pasos y en fotogramas.

        ``tiempos_fotogramas_s``: instantes de los fotogramas (p. ej. los del GIF); las
        filas se seleccionan por coincidencia de instante. Sin él, solo "todos los pasos".
        """
        bases = [("todos_los_pasos", None)]
        if tiempos_fotogramas_s is not None:
            tiempos = np.asarray(self.tiempos_s)
            marco = np.zeros(len(tiempos), dtype=bool)
            for t in np.asarray(tiempos_fotogramas_s, dtype=float):
                marco |= np.isclose(tiempos, t, rtol=0.0, atol=1e-6)
            bases.append(("fotogramas", marco))
        filas = []
        for alcance in self.alcances:
            s = self.serie(alcance)
            for base, seleccion in bases:
                candidatos = np.ones(len(s["tiempos_s"]), dtype=bool) if seleccion is None else seleccion
                for metrica, magnitud in (("dT_global_celdas", s["dt_global_c"]),
                                          ("contraste_celda_cara", s["contraste_abs_c"])):
                    validos = candidatos & np.isfinite(magnitud)
                    if not np.any(validos):
                        filas.append({"alcance": alcance, "metrica": metrica, "base": base,
                                      "valido": False})
                        continue
                    n = int(np.argmax(np.where(validos, magnitud, -np.inf)))
                    fila = {"alcance": alcance, "metrica": metrica, "base": base, "valido": True,
                            "tiempo_h": s["tiempos_s"][n] / 3600.0, "abs_c": float(magnitud[n])}
                    if metrica == "dT_global_celdas":
                        fila.update(signo_c=float(magnitud[n]), t_a_c=s["tmax"][n], t_b_c=s["tmin"][n],
                                    punto_a=self._coord(s["cel_max"][n]),    # celda más caliente
                                    punto_b=self._coord(s["cel_min"][n]),    # celda más fría
                                    cara=None)
                    else:
                        centro_cara, cara = self._coord_cara(s["con_fcel"][n], s["con_fcara"][n])
                        fila.update(signo_c=float(s["con_signo"][n]), t_a_c=s["con_tcel"][n],
                                    t_b_c=s["con_tcara"][n], punto_a=self._coord(s["con_cel"][n]),
                                    punto_b=centro_cara, cara=cara)
                    filas.append(fila)
        return filas

    # -- salidas ------------------------------------------------------------------------
    def guardar_evolucion_csv(self, ruta: Path) -> None:
        """CSV largo: una fila por instante registrado y alcance (n/a donde no hay datos)."""
        encabezado = [
            "tiempo_h", "alcance", "n_celdas_activas", "n_caras_expuestas",
            "Tmin_celdas_C", "x_Tmin_celda_m", "y_Tmin_celda_m", "z_Tmin_celda_m",
            "Tmax_celdas_C", "x_Tmax_celda_m", "y_Tmax_celda_m", "z_Tmax_celda_m", "dT_global_celdas_C",
            "Tmedia_celdas_C", "Tmin_caras_C", "x_Tmin_cara_m", "y_Tmin_cara_m", "z_Tmin_cara_m", "orientacion_Tmin_cara",
            "Tmax_caras_C", "x_Tmax_cara_m", "y_Tmax_cara_m", "z_Tmax_cara_m", "orientacion_Tmax_cara",
            "contraste_abs_celda_cara_C", "contraste_con_signo_Tcelda_menos_Tcara_C",
            "T_celda_contraste_C", "x_celda_contraste_m", "y_celda_contraste_m", "z_celda_contraste_m",
            "T_cara_contraste_C", "x_cara_contraste_m", "y_cara_contraste_m", "z_cara_contraste_m",
            "orientacion_cara_contraste"]

        def coord(plano):
            return [NA] * 3 if plano < 0 else [f"{v:.4f}" for v in self._coord(plano)]

        def coord_cara(plano, codigo):
            if plano < 0:
                return [NA] * 3, NA
            centro, cara = self._coord_cara(plano, codigo)
            return [f"{v:.4f}" for v in centro], cara

        with Path(ruta).open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow(encabezado)
            for alcance in self.alcances:
                for n, tiempo in enumerate(self.tiempos_s):
                    f = dict(zip(_CAMPOS, self._filas[alcance][n]))
                    dt = f["tmax"] - f["tmin"]
                    c_min, o_min = coord_cara(f["fmin_cel"], f["fmin_cara"])
                    c_max, o_max = coord_cara(f["fmax_cel"], f["fmax_cara"])
                    c_cara, o_con = coord_cara(f["con_fcel"], f["con_fcara"])
                    escritor.writerow([
                        f"{tiempo / 3600.0:.6f}", alcance, f["n_celdas"], f["n_caras"],
                        _fmt(f["tmin"]), *coord(f["cel_min"]), _fmt(f["tmax"]), *coord(f["cel_max"]),
                        _fmt(dt), _fmt(f["tmedia"]), _fmt(f["fmin"]), *c_min, o_min, _fmt(f["fmax"]), *c_max, o_max,
                        _fmt(abs(f["con_signo"])), _fmt(f["con_signo"]),
                        _fmt(f["con_tcel"]), *coord(f["con_cel"]), _fmt(f["con_tcara"]), *c_cara, o_con])

    def guardar_maximos_csv(self, ruta: Path, tiempos_fotogramas_s=None) -> None:
        """Tabla de máximos temporales (todos los pasos y, si se pide, solo fotogramas)."""
        with Path(ruta).open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow([
                "alcance", "metrica", "base_de_busqueda", "maximo_abs_C", "valor_con_signo_C",
                "tiempo_h", "T_punto_a_C", "x_a_m", "y_a_m", "z_a_m",
                "T_punto_b_C", "x_b_m", "y_b_m", "z_b_m", "orientacion_cara_b",
                "definicion_signo"])
            for fila in self.resumen(tiempos_fotogramas_s):
                if not fila["valido"]:
                    escritor.writerow([fila["alcance"], fila["metrica"], fila["base"], NA, NA, NA]
                                      + [NA] * 9 + [""])
                    continue
                if fila["metrica"] == "dT_global_celdas":
                    definicion = "a = celda más caliente, b = celda más fría; signo = Tmax − Tmin"
                else:
                    definicion = "a = centro de celda, b = centro de la cara; signo = Tcelda − Tcara"
                escritor.writerow([
                    fila["alcance"], fila["metrica"], fila["base"], _fmt(fila["abs_c"]),
                    f"{fila['signo_c']:+.6f}", f"{fila['tiempo_h']:.4f}",
                    _fmt(fila["t_a_c"]), *(f"{v:.4f}" for v in fila["punto_a"]),
                    _fmt(fila["t_b_c"]), *(f"{v:.4f}" for v in fila["punto_b"]),
                    fila["cara"] or NA, definicion])
