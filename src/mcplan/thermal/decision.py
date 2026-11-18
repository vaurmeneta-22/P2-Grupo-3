"""Herramientas de evaluación probabilística y criterios del dossier final.

Este módulo reutiliza ``simular_vertidas``; no incorpora una segunda ecuación de
calor. La superficie S25 se estima reconstruyendo la temperatura en la cara del
hormigón con la condición de borde actual y una interpolación lineal hasta 25 mm
hacia el interior desde el centro del voxel de borde.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from .seguimiento import CARAS


@dataclass(frozen=True)
class PlanConstruccion:
    """Calendario editable por plan; tiempos absolutos desde el inicio de obra.

    El solver actual activa cada vertida instantáneamente en ``activaciones_h``.
    ``duraciones_h`` se conserva para duración de obra y restricciones de recursos,
    pero no modela colocación parcial; esa aproximación debe declararse al reportar.
    """

    plan_id: str
    activaciones_h: dict[int, float]
    duraciones_h: dict[int, float]
    temperaturas_colocacion_c: dict[int, float] | None = None
    frentes_simultaneos_max: int | None = None

    def __post_init__(self):
        ids = set(self.activaciones_h)
        if not self.plan_id.strip() or not ids:
            raise ValueError("El plan requiere un identificador y al menos una vertida")
        if set(self.duraciones_h) != ids:
            raise ValueError("Cada vertida del plan debe tener una duración")
        if self.temperaturas_colocacion_c is not None and not set(self.temperaturas_colocacion_c) <= ids:
            raise ValueError("Temperatura asignada a una vertida que no está en el plan")
        if self.temperaturas_colocacion_c is not None and any(
                not np.isfinite(float(t)) for t in self.temperaturas_colocacion_c.values()):
            raise ValueError("Las temperaturas de colocación deben ser finitas")
        for i in ids:
            inicio, duracion = float(self.activaciones_h[i]), float(self.duraciones_h[i])
            if not np.isfinite(inicio) or inicio < 0 or not np.isfinite(duracion) or duracion < 0:
                raise ValueError(f"V{i}: activación y duración deben ser finitas y no negativas")
        if self.frentes_simultaneos_max is not None and self.frentes_simultaneos_max < 1:
            raise ValueError("frentes_simultaneos_max debe ser positivo")
        if self.frentes_simultaneos_max is not None:
            eventos = []
            for i in ids:
                inicio, fin = self.activaciones_h[i], self.activaciones_h[i] + self.duraciones_h[i]
                if fin > inicio:
                    eventos.extend(((inicio, 1), (fin, -1)))
            activos = max_activos = 0
            for _, delta in sorted(eventos, key=lambda e: (e[0], e[1])):
                activos += delta
                max_activos = max(max_activos, activos)
            if max_activos > self.frentes_simultaneos_max:
                raise ValueError("El plan excede frentes_simultaneos_max")

    def aplicar(self, vertidas: Iterable[VertidaIFC]) -> list[VertidaIFC]:
        """Devuelve copias IFC con las horas y temperaturas de este plan."""
        originales = {int(v.vertida_id): v for v in vertidas}
        if set(originales) != set(self.activaciones_h):
            raise ValueError("El plan debe definir exactamente los VertidaID del modelo IFC")
        temperaturas = self.temperaturas_colocacion_c or {}
        return [replace(originales[i], inicio_h=float(self.activaciones_h[i]),
                        temperatura_colocacion_c=float(temperaturas.get(
                            i, originales[i].temperatura_colocacion_c)))
                for i in sorted(originales)]

    def duracion_obra_h(self, t0_h: float = 0.0) -> float:
        """D = max(fin de vertida) − t0; permite que una primera colocación se retrase."""
        if not np.isfinite(t0_h):
            raise ValueError("t0_h debe ser finito")
        ultimo_fin = max(self.activaciones_h[i] + self.duraciones_h[i]
                         for i in self.activaciones_h)
        if ultimo_fin < t0_h:
            raise ValueError("El plan termina antes del t0 declarado")
        return float(ultimo_fin - t0_h)


class SeguimientoMetricasDecision:
    """Calcula M_T y M_Δ = max_t(Tmax − Tmin(S25)) en cada paso del solver.

    ``caras`` permite excluir superficies (por ejemplo, ``("z-",)`` si se decide
    no incluir la interfaz con el suelo). El valor predeterminado considera todas
    las caras exteriores tratadas por el solver, incluidas protecciones y suelo.
    """

    def __init__(self, grilla, profundidad_m: float = 0.025,
                 caras: tuple[str, ...] = CARAS, tiempo_max_s: float | None = None,
                 *, guardar_historial: bool = False):
        if not np.isfinite(profundidad_m) or profundidad_m < 0:
            raise ValueError("profundidad_m debe ser finita y no negativa")
        if not caras or any(cara not in CARAS for cara in caras):
            raise ValueError(f"caras debe ser una selección de {CARAS}")
        # ``tiempo_max_s`` limita el horizonte observado (p. ej. 14 d dentro de una
        # corrida de 28 d) sin alterar la solución: los estados posteriores se ignoran.
        self.tiempo_max_s = None if tiempo_max_s is None else float(tiempo_max_s)
        # Métricas auxiliares con nombre propio (NO sustituyen a M_delta/S25):
        # diferencia contra el centro del voxel de borde y contra la cara reconstruida.
        self.M_delta_centro_voxel_c: float | None = None
        self.M_delta_cara_c: float | None = None
        # Último estado registrado, para series temporales por paso.
        self.ultimo: dict | None = None
        self.grilla = grilla
        self.profundidad_m = float(profundidad_m)
        self.caras = tuple(caras)
        self.M_T_c: float | None = None
        self.M_delta_c: float | None = None
        self.t_M_T_s: float | None = None
        self.t_M_delta_s: float | None = None
        self.xyz_M_T_m: tuple[float, float, float] | None = None
        self.xyz_S25_min_m: tuple[float, float, float] | None = None
        self._firma_bordes = None
        self._partes = []
        self.n_muestras = 0
        self.guardar_historial = bool(guardar_historial)
        self.historial: list[dict[str, float]] = []

    def _preparar(self, solver):
        # El diccionario ``tramos_borde`` conserva su identidad cuando el solver
        # recalcula las caras tras una activación o un cambio de protecciones.
        # La máscara activa, las resistencias y las listas por cara sí cambian.
        firma = (
            id(solver.mascara_activa),
            id(solver.resistencias_borde),
            tuple(
                (cara, tuple((id(seleccion), id(condicion))
                             for seleccion, condicion in solver.tramos_borde[cara]))
                for cara in self.caras
            ),
        )
        if firma == self._firma_bordes:
            return
        partes = []
        k = solver.k.reshape(-1)
        for cara in self.caras:
            eje = "xyz".index(cara[0])
            d = float(solver.dx_m[eje])
            if self.profundidad_m > d / 2.0 + 1e-12:
                raise ValueError(
                    f"La malla en {cara} (Δx={d:g} m) no resuelve S25: se requiere Δx ≥ "
                    f"{2 * self.profundidad_m:g} m o interpolación con más capas"
                )
            for seleccion, condicion in solver.tramos_borde[cara]:
                indices = np.flatnonzero(seleccion)
                if not len(indices):
                    continue
                r_celda = d / (2.0 * k[indices])
                r_extra = solver.resistencias_borde.get(cara)
                extra = np.zeros(len(indices)) if r_extra is None else np.asarray(r_extra).reshape(-1)[indices]
                ijk = np.stack(np.unravel_index(indices, self.grilla.forma), axis=1)
                xyz = np.asarray(self.grilla.origen_m) + (ijk + 0.5) * np.asarray(self.grilla.paso_m)
                eje = "xyz".index(cara[0])
                signo_interior = -1.0 if cara.endswith("-") else 1.0
                xyz[:, eje] += signo_interior * (d / 2.0 - self.profundidad_m)
                partes.append((cara, indices, condicion, r_celda, extra, d, xyz))
        self._partes = partes
        self._firma_bordes = firma

    @staticmethod
    def _cara_c(temperatura_celda, condicion, r_celda, r_extra, tiempo_s):
        if condicion.tipo == "dirichlet":
            return np.full_like(temperatura_celda, condicion.temperatura(tiempo_s))
        if condicion.tipo == "neumann":
            return temperatura_celda + condicion.flujo_entrante(tiempo_s) * r_celda
        h = condicion.coeficiente_h_w_m2_k
        if h <= 0.0:
            return temperatura_celda.copy()
        return temperatura_celda + (condicion.temperatura(tiempo_s) - temperatura_celda) * r_celda / (
            r_celda + r_extra + 1.0 / h)

    def registrar(self, tiempo_s: float, temperatura_c, solver):
        if self.tiempo_max_s is not None and tiempo_s > self.tiempo_max_s + 1e-6:
            return
        self._preparar(solver)
        self.n_muestras += 1
        campo = np.asarray(temperatura_c)
        activa = solver.mascara_activa.reshape(-1)
        plano = campo.reshape(-1)
        if not np.any(activa):
            return
        celdas = np.flatnonzero(activa)
        i_max = int(celdas[np.argmax(plano[celdas])])
        tmax = float(plano[i_max])
        superficie_s25 = []
        coordenadas_s25 = []
        min_centro = min_cara = np.inf
        for cara, indices_original, condicion, r_original, extra_original, d, xyz_original in self._partes:
            mantener = activa[indices_original]
            if not np.any(mantener):
                continue
            indices = indices_original[mantener]
            rc, rx = r_original[mantener], extra_original[mantener]
            tc = plano[indices]
            tf = self._cara_c(tc, condicion, rc, rx, tiempo_s)
            min_centro = min(min_centro, float(np.min(tc)))
            min_cara = min(min_cara, float(np.min(tf)))
            # Distancia centro de voxel -> cara = d/2; interpolación lineal hacia adentro.
            fraccion = 0.0 if d == 0 else self.profundidad_m / (d / 2.0)
            ts25 = tf + fraccion * (tc - tf)
            superficie_s25.append(ts25)
            coordenadas_s25.append(xyz_original[mantener])
        if not superficie_s25:
            return
        todas_s25 = np.concatenate(superficie_s25)
        todas_coordenadas_s25 = np.concatenate(coordenadas_s25)
        j_min = int(np.argmin(todas_s25))
        tmin_s25 = float(todas_s25[j_min])
        delta = tmax - tmin_s25
        if self.guardar_historial:
            self.historial.append({
                "tiempo_h": float(tiempo_s) / 3600.0,
                "temperatura_maxima_c": tmax,
                "temperatura_media_c": float(np.mean(plano[celdas])),
                "temperatura_minima_voxel_c": float(np.min(plano[celdas])),
                "temperatura_minima_S25_c": tmin_s25,
                "diferencia_maxima_S25_c": delta,
            })
        self.ultimo = {"t_s": float(tiempo_s), "T_max_c": tmax, "S25_min_c": tmin_s25, "delta_c": delta}
        d_centro, d_cara = tmax - min_centro, tmax - min_cara
        if self.M_delta_centro_voxel_c is None or d_centro > self.M_delta_centro_voxel_c:
            self.M_delta_centro_voxel_c = d_centro
        if self.M_delta_cara_c is None or d_cara > self.M_delta_cara_c:
            self.M_delta_cara_c = d_cara
        xyz_max = np.unravel_index(i_max, self.grilla.forma)
        pos_max = tuple(float(self.grilla.origen_m[a] + (xyz_max[a] + 0.5) * self.grilla.paso_m[a]) for a in range(3))
        if self.M_T_c is None or tmax > self.M_T_c:
            self.M_T_c, self.t_M_T_s, self.xyz_M_T_m = tmax, float(tiempo_s), pos_max
        if self.M_delta_c is None or delta > self.M_delta_c:
            self.M_delta_c, self.t_M_delta_s = delta, float(tiempo_s)
            self.xyz_S25_min_m = tuple(float(v) for v in todas_coordenadas_s25[j_min])

    def registrar_paso(self, tiempo_s, dt_s, solver, temperatura_c, fuente_w_m3):
        """Interfaz de observador para registrar cada paso aceptado del solver."""
        self.registrar(tiempo_s, temperatura_c, solver)

    def resumen(self) -> dict:
        if self.M_T_c is None or self.M_delta_c is None:
            raise ValueError("No se registraron estados con hormigón activo y superficie S25")
        return {
            "M_T_c": self.M_T_c, "t_M_T_h": self.t_M_T_s / 3600.0,
            "xyz_M_T_m": self.xyz_M_T_m,
            "M_delta_c": self.M_delta_c, "t_M_delta_h": self.t_M_delta_s / 3600.0,
            "xyz_S25_min_m": self.xyz_S25_min_m, "profundidad_S25_m": self.profundidad_m,
            "caras_S25": list(self.caras), "n_muestras": self.n_muestras,
            "M_delta_centro_voxel_c": self.M_delta_centro_voxel_c,
            "M_delta_cara_c": self.M_delta_cara_c,
        }


@dataclass(frozen=True)
class ResultadoEscenario:
    """Una corrida plan × curva de hidratación × clima y sus métricas máximas."""

    plan_id: str
    escenario_id: str
    peso: float
    M_T_c: float
    M_delta_c: float
    t_M_T_h: float | None = None
    t_M_delta_h: float | None = None
    xyz_M_T_m: tuple[float, float, float] | None = None
    xyz_S25_min_m: tuple[float, float, float] | None = None


def agregar_probabilidades(resultados: Iterable[ResultadoEscenario], *,
                           limite_T_c: float = 70.0,
                           limite_delta_c: float = 20.0,
                           limite_P_f: float = 0.05) -> list[dict]:
    """Calcula P_T, P_Δ y la unión P_f por plan, sumando cada escenario una sola vez."""
    if (not np.isfinite(limite_T_c) or not np.isfinite(limite_delta_c)
            or not np.isfinite(limite_P_f) or not 0.0 <= limite_P_f <= 1.0):
        raise ValueError("Los límites térmicos deben ser finitos y limite_P_f debe estar entre 0 y 1")
    grupos: dict[str, list[ResultadoEscenario]] = {}
    for r in resultados:
        if not r.plan_id or not r.escenario_id:
            raise ValueError("Cada resultado necesita plan_id y escenario_id")
        if not np.isfinite(r.peso) or r.peso <= 0:
            raise ValueError("Los pesos de escenario deben ser finitos y positivos")
        if not np.isfinite(r.M_T_c) or not np.isfinite(r.M_delta_c):
            raise ValueError("Las métricas M_T y M_delta deben ser finitas")
        grupos.setdefault(r.plan_id, []).append(r)
    salida = []
    for plan_id, filas in sorted(grupos.items()):
        ids = [r.escenario_id for r in filas]
        if len(ids) != len(set(ids)):
            raise ValueError(f"El plan {plan_id} contiene escenarios duplicados")
        masa = sum(float(r.peso) for r in filas)
        if not np.isclose(masa, 1.0, rtol=0.0, atol=1e-9):
            raise ValueError(f"Los pesos del plan {plan_id} suman {masa:g}; deben sumar 1")
        falla_T = [r.M_T_c > limite_T_c for r in filas]
        falla_delta = [r.M_delta_c > limite_delta_c for r in filas]
        p_T = sum(r.peso for r, falla in zip(filas, falla_T) if falla)
        p_delta = sum(r.peso for r, falla in zip(filas, falla_delta) if falla)
        p_f = sum(r.peso for r, a, b in zip(filas, falla_T, falla_delta) if a or b)
        salida.append({"plan_id": plan_id, "n_escenarios": len(filas), "masa_pesos": masa,
                       "P_T": p_T, "P_delta": p_delta, "P_f": p_f,
                       "limite_P_f": limite_P_f, "cumple_P_f": p_f <= limite_P_f,
                       "limite_T_c": limite_T_c, "limite_delta_c": limite_delta_c})
    return salida


def guardar_resultados_escenario_csv(ruta, resultados: Iterable[ResultadoEscenario], *,
                                    limite_T_c: float = 70.0,
                                    limite_delta_c: float = 20.0) -> None:
    """Guarda trazabilidad por corrida y ubicaciones de máximos."""
    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    campos = ["plan_id", "escenario_id", "peso", "M_T_c", "M_delta_c", "t_M_T_h",
              "t_M_delta_h", "x_M_T_m", "y_M_T_m", "z_M_T_m", "x_S25_min_m",
              "y_S25_min_m", "z_S25_min_m", "I_T", "I_delta", "I_f"]
    with ruta.open("w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=campos)
        escritor.writeheader()
        for r in resultados:
            it, idelta = r.M_T_c > limite_T_c, r.M_delta_c > limite_delta_c
            coord_T = r.xyz_M_T_m or (None, None, None)
            coord_S = r.xyz_S25_min_m or (None, None, None)
            escritor.writerow({"plan_id": r.plan_id, "escenario_id": r.escenario_id, "peso": r.peso,
                               "M_T_c": r.M_T_c, "M_delta_c": r.M_delta_c,
                               "t_M_T_h": r.t_M_T_h, "t_M_delta_h": r.t_M_delta_h,
                               "x_M_T_m": coord_T[0], "y_M_T_m": coord_T[1], "z_M_T_m": coord_T[2],
                               "x_S25_min_m": coord_S[0], "y_S25_min_m": coord_S[1], "z_S25_min_m": coord_S[2],
                               "I_T": int(it), "I_delta": int(idelta), "I_f": int(it or idelta)})
