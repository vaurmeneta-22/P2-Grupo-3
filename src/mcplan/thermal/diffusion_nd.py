"""Solver explícito 2D/3D de conducción de calor por balance en cada celda.

CONVENCIONES (todas en SI; temperaturas en °C, que para diferencias es K)
--------------------------------------------------------------------------
* Grilla cartesiana regular por eje, centrada en celdas. ``T[i, j]`` (2D) o
  ``T[i, j, k]`` (3D) es la temperatura del CENTRO de la celda, en °C.
  Eje 0 = x, eje 1 = y, eje 2 = z. Forma del arreglo: ``(nx, ny[, nz])``.
* ``dx_m = (dx, dy[, dz])`` son las distancias entre centros de celdas
  vecinas, que además son los anchos de celda. La cara exterior de una celda
  de borde está a ``d/2`` de su centro.
* 2D: la celda tiene profundidad ``profundidad_m`` (1 m por defecto). Las
  energías quedan entonces en J (por la profundidad usada) y las potencias
  en W; ``Qdot`` sigue en W/m³.
* Propiedades: ``k`` [W/(m K)], ``rho`` [kg/m³], ``cp`` [J/(kg K)], escalares o
  arreglos por celda (así se representan dos materiales).

BALANCE DE UNA CELDA (esquema explícito, Euler adelantado)
----------------------------------------------------------
    rho cp V dT/dt = sum_caras Φ_cara + Qdot V       [W]

Φ_cara es el calor que ENTRA a la celda por esa cara [W]. Entre dos celdas
vecinas i y j, con resistencias en serie (media distancia en cada celda):

    R'' = d/(2 k_i) + d/(2 k_j)   [m² K/W],   Φ_i<-j = A (T_j - T_i)/R''

que equivale a usar la media armónica ``k_f = 2 k_i k_j/(k_i + k_j)``. La celda
j recibe exactamente ``-Φ`` (flujos iguales y opuestos): la energía se conserva
a precisión de máquina.

CONDICIONES DE BORDE (cada una de las 6 caras, o 4 en 2D, es independiente)
---------------------------------------------------------------------------
* Dirichlet: T prescrita EN la cara. Φ = A (T_b - T)/(d/(2k)).
* Neumann: flujo superficial ENTRANTE q_in [W/m²], Φ = q_in A. Positivo
  calienta el dominio, negativo lo enfría y cero es aislamiento. OJO: el
  solver 1D recibe el gradiente dT/dx [K/m]; aquí se recibe flujo. La
  equivalencia es q_in = -k dT/dx en la cara de x mínimo y q_in = +k dT/dx
  en la de x máximo (idem para y, z).
* Robin (convección): Φ = A (T_amb - T)/(d/(2k) + 1/h). La conducción
  media-celda y la convección están en serie. h = 0 equivale a aislamiento.

Los bordes y las esquinas salen solos: cada celda suma las caras que tiene
expuestas (una en un borde, dos en una arista 2D, tres en una esquina 3D).

ESTABILIDAD
-----------
Con Euler adelantado, T_i^{n+1} = (1 - dt ΣG_i/C_i) T_i + (positivos)·T_vecinos,
donde G [W/K] es la conductancia de cada cara (interna, Dirichlet o Robin;
Neumann no suma) y C_i = rho cp V. Los coeficientes son no negativos si

    dt <= min_i  C_i / ΣG_i

Para un interior homogéneo esto es alfa dt (1/dx² + 1/dy² [+ 1/dz²]) <= 1/2.
Las fronteras Dirichlet/Robin y las interfaces entre materiales entran en
ΣG_i, así que pueden dar un límite más estricto; el código usa el mínimo.
"""
from collections.abc import Callable
from dataclasses import dataclass, field
import math

import numpy as np

CARAS_2D = ("x-", "x+", "y-", "y+")
CARAS_3D = CARAS_2D + ("z-", "z+")
_EJES = {"x": 0, "y": 1, "z": 2}


# ---------------------------------------------------------------------------
# Condiciones de borde
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CondicionBorde:
    """Condición de una cara exterior. Usar dirichlet/neumann/robin/aislado."""

    tipo: str
    temperatura_c: float | Callable[[float], float] | None = None  # Dirichlet: T_b; Robin: T_amb
    flujo_entrante_w_m2: float | Callable[[float], float] = 0.0    # Neumann
    coeficiente_h_w_m2_k: float = 0.0                              # Robin

    def __post_init__(self):
        if self.tipo not in {"dirichlet", "neumann", "robin"}:
            raise ValueError("tipo debe ser 'dirichlet', 'neumann' o 'robin'")
        if self.tipo in {"dirichlet", "robin"} and self.temperatura_c is None:
            raise ValueError(f"{self.tipo} requiere una temperatura")
        if self.tipo == "robin":
            h = self.coeficiente_h_w_m2_k
            if not math.isfinite(h) or h < 0:
                raise ValueError("h debe ser finito y no negativo")

    def temperatura(self, tiempo_s: float) -> float:
        return _evaluar(self.temperatura_c, tiempo_s, "temperatura")

    def flujo_entrante(self, tiempo_s: float) -> float:
        return _evaluar(self.flujo_entrante_w_m2, tiempo_s, "flujo")


def dirichlet(temperatura_c) -> CondicionBorde:
    """Temperatura prescrita en la cara [°C] (constante o función de t [s])."""
    return CondicionBorde("dirichlet", temperatura_c=temperatura_c)


def neumann(flujo_entrante_w_m2=0.0) -> CondicionBorde:
    """Flujo superficial ENTRANTE al dominio [W/m²]; 0 = aislado."""
    return CondicionBorde("neumann", flujo_entrante_w_m2=flujo_entrante_w_m2)


def aislado() -> CondicionBorde:
    """Cara adiabática: Neumann con flujo cero."""
    return neumann(0.0)


def robin(coeficiente_h_w_m2_k: float, temperatura_ambiente_c) -> CondicionBorde:
    """Convección con coeficiente h [W/(m² K)] hacia un ambiente [°C]."""
    return CondicionBorde("robin", temperatura_c=temperatura_ambiente_c,
                          coeficiente_h_w_m2_k=coeficiente_h_w_m2_k)


def _evaluar(valor, tiempo_s: float, nombre: str) -> float:
    try:
        resultado = float(valor(tiempo_s) if callable(valor) else valor)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{nombre} debe ser numérico") from exc
    if not math.isfinite(resultado):
        raise ValueError(f"{nombre} debe ser finito")
    return resultado


# ---------------------------------------------------------------------------
# Resultado
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ResultadoTermico:
    """Estado final y contabilidad de energía (referencia 0 °C) de una corrida."""

    temperatura_c: np.ndarray
    tiempo_s: float
    paso_tiempo_s: float
    numero_pasos: int
    paso_maximo_estable_s: float
    energia_inicial_j: float
    energia_final_j: float
    energia_fuente_j: float
    energia_frontera_j: dict = field(default_factory=dict)  # por cara; + = entra

    @property
    def energia_frontera_total_j(self) -> float:
        return float(sum(self.energia_frontera_j.values()))

    @property
    def error_balance_j(self) -> float:
        """E_final - E_inicial - fuente - fronteras; ~0 si hay conservación."""
        return (self.energia_final_j - self.energia_inicial_j
                - self.energia_fuente_j - self.energia_frontera_total_j)


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
class DifusionCartesiana:
    """Conducción de calor 2D/3D en grilla cartesiana (volúmenes finitos)."""

    def __init__(self, *, forma, dx_m, conductividad_w_m_k, densidad_kg_m3,
                 calor_especifico_j_kg_k, fuente_w_m3=0.0, fronteras=None,
                 profundidad_m: float = 1.0):
        forma = tuple(int(n) for n in forma)
        dx_m = tuple(float(d) for d in dx_m)
        if len(forma) not in (2, 3) or len(dx_m) != len(forma):
            raise ValueError("forma y dx_m deben tener 2 o 3 elementos iguales")
        if any(n < 1 for n in forma) or any(not math.isfinite(d) or d <= 0 for d in dx_m):
            raise ValueError("n y dx deben ser positivos y finitos")
        if not math.isfinite(profundidad_m) or profundidad_m <= 0:
            raise ValueError("profundidad_m debe ser positiva y finita")
        self.forma, self.dx_m, self.ndim = forma, dx_m, len(forma)
        self.caras = CARAS_2D if self.ndim == 2 else CARAS_3D

        self.k = self._campo(conductividad_w_m_k, "k", positivo=True)
        self.rho = self._campo(densidad_kg_m3, "rho", positivo=True)
        self.cp = self._campo(calor_especifico_j_kg_k, "cp", positivo=True)
        self.fuente_w_m3 = self._campo(fuente_w_m3, "Qdot", positivo=False)

        # Volumen [m³], área de cara perpendicular a cada eje [m²], capacidad [J/K].
        espesor = (profundidad_m,) if self.ndim == 2 else ()
        self.volumen_m3 = float(np.prod(dx_m + espesor))
        self.area_cara_m2 = [float(np.prod([d for b, d in enumerate(dx_m) if b != a]
                                           + list(espesor)))
                             for a in range(self.ndim)]
        self.capacidad_j_k = self.rho * self.cp * self.volumen_m3

        self.fronteras = {cara: aislado() for cara in self.caras}
        for cara, condicion in (fronteras or {}).items():
            if cara not in self.caras:
                raise ValueError(f"cara '{cara}' inválida; usa {self.caras}")
            if not isinstance(condicion, CondicionBorde):
                raise ValueError("las fronteras deben crearse con dirichlet/neumann/robin")
            self.fronteras[cara] = condicion

        self._preparar_conductancias()

    # -- construcción -------------------------------------------------------
    def _campo(self, valor, nombre, positivo):
        arreglo = np.broadcast_to(np.asarray(valor, dtype=float), self.forma).copy()
        if not np.all(np.isfinite(arreglo)) or (positivo and np.any(arreglo <= 0)):
            raise ValueError(f"{nombre} debe ser finito"
                             + (" y positivo" if positivo else ""))
        return arreglo

    @staticmethod
    def _corte(ndim, eje, indice):
        corte = [slice(None)] * ndim
        corte[eje] = indice
        return tuple(corte)

    def _preparar_conductancias(self):
        """Conductancias G [W/K] de caras internas y de bordes Dirichlet/Robin."""
        self.G_interna = []   # G_interna[a] tiene n_a - 1 caras en el eje a
        self.G_borde = {}     # (cara) -> G [W/K] por celda de borde (Neumann: 0)
        for a in range(self.ndim):
            d, A = self.dx_m[a], self.area_cara_m2[a]
            k_bajo = self.k[self._corte(self.ndim, a, slice(None, -1))]
            k_alto = self.k[self._corte(self.ndim, a, slice(1, None))]
            # Resistencias en serie: R'' = d/2k_i + d/2k_j  (media armónica de k).
            self.G_interna.append(A / (d / (2.0 * k_bajo) + d / (2.0 * k_alto)))
            for lado, indice in (("-", 0), ("+", -1)):
                cara = "xyz"[a] + lado
                k_borde = self.k[self._corte(self.ndim, a, indice)]
                cond = self.fronteras[cara]
                resistencia = d / (2.0 * k_borde)                      # conducción
                if cond.tipo == "dirichlet":
                    G = A / resistencia
                elif cond.tipo == "robin" and cond.coeficiente_h_w_m2_k > 0:
                    G = A / (resistencia + 1.0 / cond.coeficiente_h_w_m2_k)  # + convección
                else:
                    G = np.zeros_like(k_borde)
                self.G_borde[cara] = G

    # -- estabilidad --------------------------------------------------------
    def suma_conductancias_w_k(self) -> np.ndarray:
        """ΣG de cada celda [W/K]: caras internas + Dirichlet + Robin."""
        suma = np.zeros(self.forma)
        for a in range(self.ndim):
            suma[self._corte(self.ndim, a, slice(None, -1))] += self.G_interna[a]
            suma[self._corte(self.ndim, a, slice(1, None))] += self.G_interna[a]
            for lado, indice in (("-", 0), ("+", -1)):
                suma[self._corte(self.ndim, a, indice)] += self.G_borde["xyz"[a] + lado]
        return suma

    def paso_maximo_estable_s(self) -> float:
        """dt_max = min_i C_i/ΣG_i [s] (incluye fronteras e interfaces)."""
        suma = self.suma_conductancias_w_k()
        if not np.any(suma > 0):
            return math.inf
        return float(np.min(self.capacidad_j_k[suma > 0] / suma[suma > 0]))

    def paso_estable_interior_s(self) -> float:
        """Criterio clásico del interior homogéneo: alfa dt Σ1/d² <= 1/2."""
        alfa = float(np.max(self.k / (self.rho * self.cp)))
        return 0.5 / (alfa * sum(1.0 / d**2 for d in self.dx_m))

    # -- flujos -------------------------------------------------------------
    def flujos_internos_w(self, T: np.ndarray, eje: int) -> np.ndarray:
        """Calor [W] que pasa de la celda de índice mayor a la de menor índice
        a través de cada cara interna perpendicular a ``eje`` (positivo = hacia
        -eje). Es el mismo valor, con signo opuesto, para ambas celdas."""
        alto = T[self._corte(self.ndim, eje, slice(1, None))]
        bajo = T[self._corte(self.ndim, eje, slice(None, -1))]
        return self.G_interna[eje] * (alto - bajo)

    def _flujo_borde_w(self, T, cara, tiempo_s):
        """Calor [W] que ENTRA a cada celda de la cara exterior ``cara``."""
        a = _EJES[cara[0]]
        indice = 0 if cara[1] == "-" else -1
        T_borde = T[self._corte(self.ndim, a, indice)]
        cond = self.fronteras[cara]
        if cond.tipo == "neumann":                       # Φ = q_in A
            return np.full(T_borde.shape,
                           cond.flujo_entrante(tiempo_s) * self.area_cara_m2[a])
        # Dirichlet: T_b ; Robin: T_amb. Ambos: Φ = G (T_ref - T)
        return self.G_borde[cara] * (cond.temperatura(tiempo_s) - T_borde)

    def _derivada_energia_w(self, T, tiempo_s):
        """Potencia neta que entra a cada celda [W] y por cara exterior."""
        potencia = self.fuente_w_m3 * self.volumen_m3            # Qdot V
        for a in range(self.ndim):
            f = self.flujos_internos_w(T, a)
            potencia[self._corte(self.ndim, a, slice(None, -1))] += f   # entra a la celda baja
            potencia[self._corte(self.ndim, a, slice(1, None))] -= f    # sale de la alta
        entrada_frontera = {}
        for cara in self.caras:
            flujo = self._flujo_borde_w(T, cara, tiempo_s)
            potencia[self._corte(self.ndim, _EJES[cara[0]],
                                 0 if cara[1] == "-" else -1)] += flujo
            entrada_frontera[cara] = float(flujo.sum())
        return potencia, entrada_frontera

    def energia_j(self, T: np.ndarray) -> float:
        """Energía interna respecto de 0 °C: Σ rho cp V T [J]."""
        return float(np.sum(self.capacidad_j_k * T))

    # -- integración --------------------------------------------------------
    def simular(self, *, temperatura_inicial_c, tiempo_final_s: float,
                paso_tiempo_max_s: float) -> ResultadoTermico:
        """Avanza hasta ``tiempo_final_s`` con pasos iguales <= paso_tiempo_max_s.

        Rechaza (ValueError) pasos mayores que el límite de estabilidad.
        """
        if not all(math.isfinite(v) and v > 0 for v in (tiempo_final_s, paso_tiempo_max_s)):
            raise ValueError("tiempo final y paso deben ser positivos y finitos")
        T = self._campo(temperatura_inicial_c, "T inicial", positivo=False)
        numero_pasos = math.ceil(tiempo_final_s / paso_tiempo_max_s)
        dt = tiempo_final_s / numero_pasos
        dt_max = self.paso_maximo_estable_s()
        if dt > dt_max * (1.0 + 1e-12):
            raise ValueError(
                f"Paso inestable: dt={dt:.6g} s supera el límite {dt_max:.6g} s "
                f"(interior homogéneo: {self.paso_estable_interior_s():.6g} s)")

        energia_inicial = self.energia_j(T)
        energia_fuente = float(np.sum(self.fuente_w_m3 * self.volumen_m3)) * tiempo_final_s
        energia_frontera = {cara: 0.0 for cara in self.caras}
        for n in range(numero_pasos):
            potencia, entrada = self._derivada_energia_w(T, n * dt)
            T = T + dt * potencia / self.capacidad_j_k
            for cara, w in entrada.items():
                energia_frontera[cara] += w * dt
        return ResultadoTermico(
            temperatura_c=T, tiempo_s=tiempo_final_s, paso_tiempo_s=dt,
            numero_pasos=numero_pasos, paso_maximo_estable_s=dt_max,
            energia_inicial_j=energia_inicial, energia_final_j=self.energia_j(T),
            energia_fuente_j=energia_fuente, energia_frontera_j=energia_frontera)


def interpolar_en_punto(temperatura_c, dx_m, punto_m) -> float:
    """Temperatura en un punto físico [m] por interpolación multilineal.

    Los centros de celda están en (i + 1/2)·d; fuera de los centros extremos
    se mantiene el valor de la celda de borde. Sirve para comparar mallas
    distintas en el MISMO punto (T[n//2] cambia de posición al refinar).
    """
    T = np.asarray(temperatura_c, dtype=float)
    if T.ndim != len(dx_m) or len(punto_m) != T.ndim:
        raise ValueError("dimensiones inconsistentes")
    for eje, (d, x) in enumerate(zip(dx_m, punto_m)):
        n = T.shape[0]
        s = min(max(x / d - 0.5, 0.0), n - 1.0)
        i = min(int(math.floor(s)), max(n - 2, 0))
        w = s - i if n > 1 else 0.0
        T = (1.0 - w) * T[i] + w * T[min(i + 1, n - 1)]
    return float(T)
