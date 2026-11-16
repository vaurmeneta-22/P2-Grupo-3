"""Registro de la generación por hidratación y de los flujos hacia el entorno, paso a paso.

Observador de ``simular_vertidas`` (parámetro ``seguimiento``). Lee, en cada paso de integración,
los MISMOS valores que usa el solver (no recalcula nada ni guarda campos completos):

* ``Qdot`` [W/m³]: el campo de generación del paso (``solver._fuente_en_tiempo``), por vertida,
  y la potencia ``Σ Qdot·V_voxel`` [W], la misma que acumula el solver en ``energia_fuente_j``.
* Flujos por cara exterior: la potencia [W] que entra por cada celda
  (``solver.flujos_ultimo_paso_w``), la misma que el solver suma en ``energia_frontera_j``.

Convenio de signo: **positivo = entra al hormigón, negativo = sale hacia el entorno**.

Unidades y muestreo:
* Las muestras corresponden al INSTANTE INICIAL de cada paso (el solver es explícito: usa el campo
  y los valores de ese instante para avanzar el paso).
* Flujo por unidad de área ``q`` [W/m²] = potencia de las celdas seleccionadas / (n_celdas · área de cara).
* Potencia integrada [W] = suma de la potencia de esas celdas (depende de cuántas celdas expuestas haya,
  por eso no es comparable entre vertidas de distinto tamaño; ``q`` sí compara intensidad).
* El suelo se registra aparte: celdas de la base física con cara ``z-`` expuesta.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def mascara_base(grilla) -> np.ndarray:
    """Celdas de la base física: la capa más baja del hormigón (misma regla que ``ModeloBordes``)."""
    solido = grilla.vertida_id > 0
    capas = np.flatnonzero(solido.any(axis=(0, 1)))
    k0 = int(capas[0]) if capas.size else 0
    return solido & (np.arange(grilla.forma[2])[None, None, :] == k0)


class RegistroEnergetico:
    """Observador de pasos: Qdot por vertida, potencia generada y flujos por caras y suelo."""

    def __init__(self, grilla, *, caras=(), celdas_suelo=None):
        """``caras``: pares ``(vertida_id, cara)`` representativos, p. ej. ``[(1, "x+"), (5, "x+")]``.

        ``celdas_suelo``: máscara booleana de las celdas de la base física (``ModeloBordes.mascara_base``).
        """
        self.grilla = grilla
        self.ids = sorted(int(i) for i in np.unique(grilla.vertida_id) if i > 0)
        self.caras = [(int(v), str(c)) for v, c in caras]
        self._celdas = {i: grilla.vertida_id == i for i in self.ids}
        self.celdas_suelo = None if celdas_suelo is None else np.asarray(celdas_suelo, dtype=bool)
        self.t_s: list[float] = []
        self.dt_s: list[float] = []
        self.potencia_total_w: list[float] = []
        self.qdot_w_m3 = {i: [] for i in self.ids}
        self.potencia_vertida_w = {i: [] for i in self.ids}
        self.frontera_total_w: list[float] = []          # Σ de todas las caras (comparable con el solver)
        self.suelo_w: list[float] = []
        self.suelo_w_m2: list[float] = []
        self.suelo_n: list[int] = []
        self.cara_w = {k: [] for k in self.caras}
        self.cara_w_m2 = {k: [] for k in self.caras}
        self.cara_n = {k: [] for k in self.caras}
        self.cara_resistencia = {k: [] for k in self.caras}   # R'' de protección vigente [m² K/W]
        self.energia_fuente_j = 0.0
        self.energia_frontera_j = 0.0

    # interfaz de observador de la simulación (``registrar`` no necesita nada: todo viene por paso)
    def registrar(self, tiempo_s, temperatura_c, solver) -> None:
        return None

    def registrar_paso(self, t_s, dt_s, solver, temperatura_c, fuente_w_m3) -> None:
        volumen = solver.volumen_m3
        self.t_s.append(float(t_s))
        self.dt_s.append(float(dt_s))
        potencia = float(np.sum(fuente_w_m3 * volumen))
        self.potencia_total_w.append(potencia)
        self.energia_fuente_j += potencia * dt_s
        activa = solver.mascara_activa
        for i in self.ids:
            celdas = self._celdas[i]
            en_vertida = celdas & activa
            self.qdot_w_m3[i].append(float(np.mean(fuente_w_m3[en_vertida])) if np.any(en_vertida) else math.nan)
            self.potencia_vertida_w[i].append(float(np.sum(fuente_w_m3[celdas]) * volumen)
                                              if np.any(en_vertida) else math.nan)
        flujos = solver.flujos_ultimo_paso_w
        total = sum(float(f.sum()) for f in flujos.values())
        self.frontera_total_w.append(total)
        self.energia_frontera_j += total * dt_s
        if self.celdas_suelo is not None:
            celdas = self.celdas_suelo & solver.superficie_expuesta["z-"]
            n = int(celdas.sum())
            w = float(flujos["z-"][celdas].sum()) if n else math.nan
            self.suelo_n.append(n)
            self.suelo_w.append(w)
            self.suelo_w_m2.append(w / (n * solver.area_cara_m2[2]) if n else math.nan)
        for clave in self.caras:
            vertida, cara = clave
            celdas = solver.superficie_expuesta[cara] & self._celdas[vertida]
            n = int(celdas.sum())
            area = solver.area_cara_m2["xyz".index(cara[0])]
            w = float(flujos[cara][celdas].sum()) if n else math.nan
            r = solver.resistencias_borde.get(cara)
            self.cara_n[clave].append(n)
            self.cara_w[clave].append(w)
            self.cara_w_m2[clave].append(w / (n * area) if n else math.nan)
            self.cara_resistencia[clave].append(float(np.mean(np.asarray(r)[celdas])) if (n and r is not None) else
                                                (0.0 if n else math.nan))

    # -- consulta ------------------------------------------------------------------
    @property
    def horas(self) -> np.ndarray:
        return np.asarray(self.t_s) / 3600.0

    def potencia_total_kw(self) -> np.ndarray:
        return np.asarray(self.potencia_total_w) / 1000.0

    def pico(self, serie, signo=1.0) -> tuple[float, float] | None:
        """(hora [h], valor) del máximo de ``signo * serie`` ignorando NaN; ``None`` si no hay datos."""
        y = np.asarray(serie, dtype=float)
        if not np.any(np.isfinite(y)):
            return None
        n = int(np.nanargmax(signo * y))
        return float(self.horas[n]), float(y[n])

    # -- CSV ------------------------------------------------------------------------
    def guardar_hidratacion_csv(self, ruta) -> None:
        with Path(ruta).open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["tiempo_h", "potencia_total_W", "potencia_total_kW"]
                       + [f"Qdot_V{i}_W_m3" for i in self.ids] + [f"potencia_V{i}_W" for i in self.ids])
            for n, t in enumerate(self.horas):
                w.writerow([f"{t:.6f}", f"{self.potencia_total_w[n]:.6f}", f"{self.potencia_total_w[n] / 1000:.9f}"]
                           + [_f(self.qdot_w_m3[i][n]) for i in self.ids]
                           + [_f(self.potencia_vertida_w[i][n]) for i in self.ids])

    def guardar_flujos_csv(self, ruta) -> None:
        """Una fila por paso. Convenio: + entra al hormigón, − sale. W/m² = flujo por área; W = potencia."""
        with Path(ruta).open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            cabecera = ["tiempo_h", "frontera_total_W"]
            if self.celdas_suelo is not None:
                cabecera += ["suelo_W_m2", "suelo_W", "suelo_n_celdas"]
            for v, c in self.caras:
                cabecera += [f"V{v}_{c}_W_m2", f"V{v}_{c}_W", f"V{v}_{c}_n_celdas", f"V{v}_{c}_R_proteccion_m2K_W"]
            w.writerow(cabecera)
            for n, t in enumerate(self.horas):
                fila = [f"{t:.6f}", f"{self.frontera_total_w[n]:.6f}"]
                if self.celdas_suelo is not None:
                    fila += [_f(self.suelo_w_m2[n]), _f(self.suelo_w[n]), self.suelo_n[n]]
                for clave in self.caras:
                    fila += [_f(self.cara_w_m2[clave][n]), _f(self.cara_w[clave][n]), self.cara_n[clave][n],
                             _f(self.cara_resistencia[clave][n])]
                w.writerow(fila)


def _f(valor) -> str:
    return "n/a" if valor is None or not math.isfinite(valor) else f"{valor:.8g}"
