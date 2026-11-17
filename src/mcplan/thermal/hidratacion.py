"""Modelos de fuente volumétrica dependientes de la edad de una vertida.

Incluye una curva adiabática exponencial y el pulso cosenoidal usado por el
pipeline IFC 3D. Sus parámetros son configurables y no son una calibración
del hormigón real del proyecto.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ModeloHidratacionExponencial:
    """Elevación adiabática asintótica y tiempo característico.

    ``T_ad(age) = T_colocación + delta_T * (1-exp(-age/tau))``.
    La fuente se deriva como ``Qdot = rho*cp*dT_ad/dt`` [W/m³].
    """

    incremento_adiabatico_c: float
    tiempo_caracteristico_s: float

    def __post_init__(self):
        if not math.isfinite(self.incremento_adiabatico_c) or self.incremento_adiabatico_c < 0:
            raise ValueError("incremento_adiabatico_c debe ser finito y no negativo")
        if not math.isfinite(self.tiempo_caracteristico_s) or self.tiempo_caracteristico_s <= 0:
            raise ValueError("tiempo_caracteristico_s debe ser positivo y finito")

    def temperatura_adiabatica_c(self, temperatura_colocacion_c: float,
                                 edad_s) -> np.ndarray | float:
        edad = np.maximum(np.asarray(edad_s, dtype=float), 0.0)
        temperatura = temperatura_colocacion_c + self.incremento_adiabatico_c * (
            1.0 - np.exp(-edad / self.tiempo_caracteristico_s))
        return float(temperatura) if temperatura.ndim == 0 else temperatura

    def generacion_w_m3(self, edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Devuelve Qdot [W/m³], nula antes de colocar y positiva después."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(rho <= 0) or np.any(cp <= 0):
            raise ValueError("rho y cp deben ser positivos")
        fuente = rho * cp * self.incremento_adiabatico_c / self.tiempo_caracteristico_s
        return np.where(edad >= 0.0, fuente * np.exp(-np.maximum(edad, 0.0)
                                                      / self.tiempo_caracteristico_s), 0.0)


@dataclass(frozen=True)
class ModeloPulsoCosenoidal:
    """Fuente volumétrica cosenoidal definida por la edad de cada vertida.

    ``Q(age) = Qmax/2 * (1-cos(2*pi*age/duracion))`` para
    ``0 <= age <= duracion`` y cero fuera de ese intervalo. La integral
    adiabática del pulso es ``Qmax*duracion/(2*rho*cp)`` [K].
    """

    qmax_w_m3: float
    duracion_s: float

    def __post_init__(self):
        if not math.isfinite(self.qmax_w_m3) or self.qmax_w_m3 < 0:
            raise ValueError("qmax_w_m3 debe ser finito y no negativo")
        if not math.isfinite(self.duracion_s) or self.duracion_s <= 0:
            raise ValueError("duracion_s debe ser positiva y finita")

    def generacion_w_m3(self, edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Devuelve Q(age) en W/m³; cada elemento puede tener edad distinta."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
            raise ValueError("rho debe ser finita y positiva")
        if np.any(~np.isfinite(cp)) or np.any(cp <= 0):
            raise ValueError("cp debe ser finito y positivo")
        activa = (edad >= 0.0) & (edad <= self.duracion_s)
        edad_segura = np.clip(edad, 0.0, self.duracion_s)
        q = 0.5 * self.qmax_w_m3 * (
            1.0 - np.cos(2.0 * math.pi * edad_segura / self.duracion_s)
        )
        return np.where(activa, q, 0.0)

    def temperatura_adiabatica_c(self, temperatura_colocacion_c: float,
                                  edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Temperatura acumulada ideal sin pérdidas para una vertida [°C]."""
        edad = np.asarray(edad_s, dtype=float)
        rho = np.asarray(densidad_kg_m3, dtype=float)
        cp = np.asarray(calor_especifico_j_kg_k, dtype=float)
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
            raise ValueError("rho debe ser finita y positiva")
        if np.any(~np.isfinite(cp)) or np.any(cp <= 0):
            raise ValueError("cp debe ser finito y positivo")
        edad_limitada = np.clip(edad, 0.0, self.duracion_s)
        energia_por_volumen = 0.5 * self.qmax_w_m3 * (
            edad_limitada - self.duracion_s / (2.0 * math.pi)
            * np.sin(2.0 * math.pi * edad_limitada / self.duracion_s)
        )
        temperatura = temperatura_colocacion_c + energia_por_volumen / (rho * cp)
        return float(temperatura) if temperatura.ndim == 0 else temperatura


@dataclass(frozen=True)
class ModeloCurvaAdabatica:
    """Fuente volumétrica obtenida de una curva adiabática medida/tabulada.

    ``tiempos_s`` son edades desde la colocación y ``incrementos_c`` son
    elevaciones sobre la temperatura de colocación. Entre puntos, la curva es
    lineal; por tanto, su derivada (y ``Qdot``) es constante por tramo. Fuera
    del intervalo medido, la tasa es cero. No extrapola ni crea una cola ficticia.
    """

    tiempos_s: tuple[float, ...]
    incrementos_c: tuple[float, ...]

    def __post_init__(self):
        t = np.asarray(self.tiempos_s, dtype=float)
        d = np.asarray(self.incrementos_c, dtype=float)
        if t.ndim != 1 or d.ndim != 1 or len(t) != len(d) or len(t) < 2:
            raise ValueError("La curva adiabática requiere al menos dos pares tiempo-temperatura")
        if not np.all(np.isfinite(t)) or not np.all(np.isfinite(d)):
            raise ValueError("La curva adiabática debe contener solo valores finitos")
        if t[0] < 0 or np.any(np.diff(t) <= 0):
            raise ValueError("Las edades deben empezar en 0 o después y ser estrictamente crecientes")
        if np.any(np.diff(d) < -1e-10):
            raise ValueError("La elevación adiabática debe ser monótona no decreciente")

    @classmethod
    def desde_csv(cls, ruta, columna_tiempo="edad_h", columna_incremento="incremento_adiabatico_c",
                  ajuste_monotono="ninguno"):
        """Carga edad [h] y elevación [°C]; puede aplicar isotonic regression explícita.

        ``ajuste_monotono="isotonico"`` corrige ruido de medición a la curva
        no decreciente de mínimos cuadrados. El CSV original no se modifica.
        """
        import csv
        from pathlib import Path

        with Path(ruta).open(newline="", encoding="utf-8-sig") as archivo:
            lector = csv.DictReader(archivo)
            if not lector.fieldnames or columna_tiempo not in lector.fieldnames or columna_incremento not in lector.fieldnames:
                raise ValueError(f"CSV debe incluir columnas {columna_tiempo!r} y {columna_incremento!r}")
            filas = list(lector)
        try:
            tiempos = tuple(float(f[columna_tiempo]) * 3600.0 for f in filas)
            incrementos = tuple(float(f[columna_incremento]) for f in filas)
        except (TypeError, ValueError) as exc:
            raise ValueError("La curva CSV contiene valores no numéricos") from exc
        ajuste = str(ajuste_monotono).lower()
        if ajuste not in ("ninguno", "isotonico"):
            raise ValueError("ajuste_monotono debe ser 'ninguno' o 'isotonico'")
        if ajuste == "isotonico":
            incrementos = cls._ajuste_isotonico(incrementos)
        return cls(tiempos, incrementos)

    @staticmethod
    def _ajuste_isotonico(valores):
        """Proyección PAVA de valores a una secuencia monótona no decreciente."""
        bloques = []
        for valor in map(float, valores):
            bloques.append([valor, 1])
            while len(bloques) >= 2 and bloques[-2][0] > bloques[-1][0]:
                v1, n1 = bloques[-2]
                v2, n2 = bloques[-1]
                bloques[-2:] = [[(v1 * n1 + v2 * n2) / (n1 + n2), n1 + n2]]
        return tuple(valor for valor, n in bloques for _ in range(n))

    def generacion_w_m3(self, edad_s, densidad_kg_m3, calor_especifico_j_kg_k):
        """Evalúa ``rho·cp·d(ΔT_ad)/dt`` en W/m³ para edades escalares o vectoriales."""
        edad, rho, cp = np.broadcast_arrays(
            np.asarray(edad_s, dtype=float), np.asarray(densidad_kg_m3, dtype=float),
            np.asarray(calor_especifico_j_kg_k, dtype=float))
        if np.any(~np.isfinite(rho)) or np.any(rho <= 0) or np.any(~np.isfinite(cp)) or np.any(cp <= 0):
            raise ValueError("rho y cp deben ser finitos y positivos")
        tiempos = np.asarray(self.tiempos_s, dtype=float)
        incrementos = np.asarray(self.incrementos_c, dtype=float)
        pendientes = np.diff(incrementos) / np.diff(tiempos)
        tramo = np.searchsorted(tiempos, edad, side="right") - 1
        dentro = (edad >= tiempos[0]) & (edad < tiempos[-1])
        pendiente = np.zeros(edad.shape, dtype=float)
        indices = np.clip(tramo, 0, len(pendientes) - 1)
        pendiente[dentro] = pendientes[indices[dentro]]
        salida = rho * cp * pendiente
        return float(salida) if salida.ndim == 0 else salida

    def temperatura_adiabatica_c(self, temperatura_colocacion_c: float, edad_s):
        """Interpola la curva adiabática; después del último dato mantiene el último nivel."""
        edad = np.asarray(edad_s, dtype=float)
        valor = np.interp(np.maximum(edad, 0.0), self.tiempos_s, self.incrementos_c,
                          left=0.0, right=self.incrementos_c[-1]) + temperatura_colocacion_c
        return float(valor) if valor.ndim == 0 else valor


@dataclass(frozen=True)
class ModeloFuenteTabulada:
    """Fuente ``Qdot(edad)`` tabulada directamente en W/m³.

    Interpola linealmente entre muestras. Con ``cola_cero=True``, fija la fuente
    en cero después del último dato válido; esa extrapolación queda explícita y
    se usa para curvas con ventana de registro finita.
    """

    tiempos_s: tuple[float, ...]
    qdot_w_m3: tuple[float, ...]
    cola_cero: bool = True

    def __post_init__(self):
        t = np.asarray(self.tiempos_s, dtype=float)
        q = np.asarray(self.qdot_w_m3, dtype=float)
        if t.ndim != 1 or q.ndim != 1 or len(t) != len(q) or len(t) < 2:
            raise ValueError("La fuente tabulada requiere al menos dos pares tiempo-Qdot")
        if not np.all(np.isfinite(t)) or not np.all(np.isfinite(q)):
            raise ValueError("La fuente tabulada debe ser finita; recorta NaN finales explícitamente")
        if t[0] < 0 or np.any(np.diff(t) <= 0):
            raise ValueError("Las edades deben empezar en 0 o después y ser estrictamente crecientes")
        if np.any(q < 0):
            raise ValueError("Qdot tabulado debe ser no negativo")

    @classmethod
    def desde_source_terms(cls, ruta, columna: str, *, cola_cero: bool = True):
        """Lee una columna de `source_terms.txt`; admite NaN solo al final de esa curva."""
        from pathlib import Path

        lineas = Path(ruta).read_text(encoding="utf-8-sig").splitlines()
        cabecera = None
        filas = []
        for linea in lineas:
            texto = linea.strip()
            if not texto:
                continue
            if texto.startswith("#"):
                posible = texto.lstrip("# ").split()
                if posible and posible[0] == "time_h":
                    cabecera = posible
                continue
            if cabecera is None:
                continue
            valores = texto.split()
            if len(valores) != len(cabecera):
                raise ValueError(f"Fila mal formada en source_terms: {texto}")
            filas.append(valores)
        if cabecera is None or columna not in cabecera:
            raise ValueError(f"source_terms debe tener una columna {columna!r}")
        if len(filas) < 2:
            raise ValueError("source_terms no contiene suficientes muestras")
        indice = cabecera.index(columna)
        tiempos_h, qdot = [], []
        faltante_visto = False
        for fila in filas:
            try:
                t_h = float(fila[0])
                q = float(fila[indice])
            except ValueError as exc:
                raise ValueError(f"Valor no numérico en source_terms: {fila}") from exc
            if not np.isfinite(t_h):
                raise ValueError("Las edades de source_terms deben ser finitas")
            if np.isnan(q):
                faltante_visto = True
                continue
            if faltante_visto:
                raise ValueError(f"La columna {columna} tiene datos después de un NaN intermedio")
            tiempos_h.append(t_h)
            qdot.append(q)
        if len(tiempos_h) < 2:
            raise ValueError(f"La curva {columna} tiene menos de dos puntos observados")
        return cls(tuple(t * 3600.0 for t in tiempos_h), tuple(qdot), cola_cero)

    def generacion_w_m3(self, edad_s, densidad_kg_m3=None, calor_especifico_j_kg_k=None):
        """Evalúa Qdot [W/m³]; acepta los argumentos rho/cp del protocolo de hidratación."""
        edad = np.asarray(edad_s, dtype=float)
        if not self.cola_cero and np.any(edad > self.tiempos_s[-1]):
            raise ValueError("La edad excede la cobertura de la curva tabulada")
        q = np.interp(edad, self.tiempos_s, self.qdot_w_m3,
                      left=0.0, right=0.0 if self.cola_cero else self.qdot_w_m3[-1])
        q = np.where(edad >= 0.0, q, 0.0)
        return float(q) if q.ndim == 0 else q
