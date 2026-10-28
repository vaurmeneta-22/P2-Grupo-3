"""Pruebas locales del ejercicio de fuente variable (solo usan solver_ejercicio.py de esta carpeta).

Ejecutar:  python -m unittest discover -s ejercicios_locales/fuente_variable_2d -v
(o pulsar ▶️ en este archivo). No forman parte de la bateria de la entrega (tests/).
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import math
import os
import sys
import unittest

os.environ.setdefault("MPLBACKEND", "Agg")      # el script importa matplotlib
import numpy as np

CARPETA = Path(__file__).resolve().parent


def _cargar(nombre, archivo):
    spec = spec_from_file_location(nombre, CARPETA / archivo)
    modulo = module_from_spec(spec)
    sys.modules[nombre] = modulo
    spec.loader.exec_module(modulo)
    return modulo


solver = _cargar("solver_ejercicio", "solver_ejercicio.py")

K, RHO, CP = 2.0, 2400.0, 900.0
RHOCP = RHO * CP


def q_pulso(t, qmax, tau):
    """Q(t) = Qmax/2 [1 - cos(2 pi t/tau)] para 0 <= t <= tau; 0 fuera."""
    return 0.5 * qmax * (1.0 - math.cos(2.0 * math.pi * t / tau)) if 0.0 <= t <= tau else 0.0


def sistema(fuente, fronteras=None, forma=(6, 4), dx=(0.1, 0.1)):
    return solver.DifusionCartesiana(
        forma=forma, dx_m=dx, conductividad_w_m_k=K, densidad_kg_m3=RHO,
        calor_especifico_j_kg_k=CP, fuente_w_m3=fuente, fronteras=fronteras)


def aislados():
    return {c: solver.aislado() for c in ("x-", "x+", "y-", "y+")}


def dirichlet20():
    return {c: solver.dirichlet(20.0) for c in ("x-", "x+", "y-", "y+")}


class TestFuenteConstanteComoFuncion(unittest.TestCase):
    def test_constante_escalar_igual_a_funcion_constante(self):
        for q in (0.0, 250.0):
            a = sistema(q, aislados())
            b = sistema(lambda t, q=q: q, aislados())
            kw = dict(temperatura_inicial_c=20.0, tiempo_final_s=7200.0,
                      paso_tiempo_max_s=a.paso_maximo_estable_s())
            ra, rb = a.simular(**kw), b.simular(**kw)
            np.testing.assert_array_equal(ra.temperatura_c, rb.temperatura_c)
            self.assertAlmostEqual(ra.energia_fuente_j, rb.energia_fuente_j, delta=1e-9 * (1 + abs(ra.energia_fuente_j)))

    def test_calentamiento_uniforme_aislado(self):
        q, t = 500.0, 6 * 3600.0
        s = sistema(lambda _t: q, aislados())
        r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=t, paso_tiempo_max_s=s.paso_maximo_estable_s())
        np.testing.assert_allclose(r.temperatura_c, 20.0 + q * t / RHOCP, rtol=1e-12)


class TestLinealidad(unittest.TestCase):
    def test_escalamiento_lineal_con_qmax(self):
        tau, T = 7200.0, 7200.0
        centro = [(0.3, 0.2)]
        res = {}
        for qmax in (50.0, 150.0, 400.0):
            s = sistema(lambda t, q=qmax: q_pulso(t, q, tau), dirichlet20())
            res[qmax] = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=T,
                                  paso_tiempo_max_s=0.8 * s.paso_maximo_estable_s(), puntos_registro_m=centro)
        base = res[50.0]
        for qmax in (150.0, 400.0):
            f = qmax / 50.0
            np.testing.assert_allclose(res[qmax].temperatura_c - 20.0, f * (base.temperatura_c - 20.0), rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(res[qmax].historial_puntos_c[:, 0] - 20.0,
                                       f * (base.historial_puntos_c[:, 0] - 20.0), rtol=1e-10, atol=1e-12)

    def test_qmax_por_linealidad_alcanza_el_limite(self):
        tau, T, limite = 7200.0, 7200.0, 70.0
        def maximo(qmax):
            s = sistema(lambda t: q_pulso(t, qmax, tau), dirichlet20())
            r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=T,
                          paso_tiempo_max_s=0.8 * s.paso_maximo_estable_s(), puntos_registro_m=[(0.3, 0.2)])
            return r.historial_puntos_c[:, 0].max()
        qref = 80.0
        qlim = qref * (limite - 20.0) / (maximo(qref) - 20.0)
        self.assertAlmostEqual(maximo(qlim), limite, delta=1e-9)
        self.assertLess(maximo(0.999 * qlim), limite)
        self.assertGreater(maximo(1.001 * qlim), limite)


class TestBalanceConFuenteTemporal(unittest.TestCase):
    def test_residuo_y_energia_de_la_fuente(self):
        tau, T, qmax = 7200.0, 7200.0, 300.0
        s = sistema(lambda t: q_pulso(t, qmax, tau), dirichlet20())
        n = 200
        r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=T, paso_tiempo_max_s=T / n * (1 + 1e-9))
        self.assertEqual(r.numero_pasos, n)
        self.assertLess(abs(r.error_balance_j), 1e-9 * abs(r.energia_final_j))
        # integral exacta de Q V sobre un periodo completo = Qmax V tau / 2 (la suma por la izquierda
        # de una funcion periodica suave es exacta a precision de maquina cuando tau = n dt)
        vol_total = s.volumen_m3 * 6 * 4
        self.assertAlmostEqual(r.energia_fuente_j, qmax * vol_total * tau / 2.0, delta=1e-9 * qmax * vol_total * tau)
        self.assertGreater(r.energia_frontera_total_j, -1.0 * r.energia_fuente_j)   # algo sale por los bordes

    def test_historial_en_todos_los_pasos(self):
        s = sistema(lambda t: 100.0, aislados())
        r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=3600.0, paso_tiempo_max_s=100.0,
                      puntos_registro_m=[(0.3, 0.2)])
        self.assertEqual(len(r.historial_tiempos_s), r.numero_pasos + 1)
        self.assertEqual(r.historial_puntos_c.shape, (r.numero_pasos + 1, 1))
        self.assertAlmostEqual(r.historial_tiempos_s[-1], 3600.0)
        self.assertEqual(r.historial_puntos_c[0, 0], 20.0)
        self.assertAlmostEqual(r.historial_puntos_c[-1, 0], 20.0 + 100.0 * 3600.0 / RHOCP, places=9)

    def test_interpolacion_precalculada_coincide_con_interpolar_en_punto(self):
        rng = np.random.default_rng(0)
        s = sistema(0.0, aislados(), forma=(8, 6), dx=(0.1, 0.2))
        T = rng.uniform(0, 50, (8, 6))
        for p in ((0.4, 0.6), (0.0, 0.0), (0.8, 1.2), (0.03, 1.19), (0.55, 0.3)):
            w = s._pesos_interpolacion(p)
            self.assertAlmostEqual(s._evaluar_pesos(T, w), solver.interpolar_en_punto(T, s.dx_m, p), places=12)
        with self.assertRaises(ValueError):
            s._pesos_interpolacion((5.0, 0.1))


class TestPulsoAdiabaticoAnalitico(unittest.TestCase):
    """DeltaT(t) = Qmax/(2 rho cp) [t - tau/(2pi) sin(2 pi t/tau)] (t <= tau); Qmax tau/(2 rho cp) despues.

    Tolerancia (error de integracion temporal): el solver acumula por la IZQUIERDA,
    sum_{k<n} Q(t_k) dt. Con Euler-Maclaurin, para Q(0) = 0 y Q'(0) = 0:
        DeltaT_num(t_n) - DeltaT_exacta(t_n) = -(dt/2) Q(t_n)/(rho cp) + O(dt^2 Q'/(rho cp)).
    Cota usada: (dt/2) Qmax/(rho cp) + (dt²/12)(pi Qmax/tau)/(rho cp)  (max |Q'| = pi Qmax/tau).
    Para t >= tau y tau = N dt la suma por la izquierda de un periodo completo es exacta.
    """

    def exacto(self, t, qmax, tau):
        t = np.minimum(t, tau)
        return qmax / (2.0 * RHOCP) * (t - tau / (2 * math.pi) * np.sin(2 * math.pi * t / tau))

    def correr(self, n_por_tau, horizonte_en_tau=1.5, qmax=500.0, tau=7200.0):
        s = sistema(lambda t: q_pulso(t, qmax, tau), aislados(), forma=(4, 4))
        dt = tau / n_por_tau
        T = horizonte_en_tau * tau
        r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=T, paso_tiempo_max_s=dt * (1 + 1e-9),
                      puntos_registro_m=[(0.2, 0.2)])
        return r, dt, qmax, tau

    def test_serie_completa_dentro_de_la_tolerancia(self):
        for n_por_tau in (100, 200, 400):
            r, dt, qmax, tau = self.correr(n_por_tau)
            t = r.historial_tiempos_s
            err = (r.historial_puntos_c[:, 0] - 20.0) - self.exacto(t, qmax, tau)
            tol = 0.5 * dt * qmax / RHOCP + dt**2 / 12.0 * (math.pi * qmax / tau) / RHOCP
            self.assertLessEqual(np.max(np.abs(err)), tol, (n_por_tau, np.max(np.abs(err)), tol))
            # el termino principal predice el error dentro del orden dt²
            pred = -0.5 * dt * np.array([q_pulso(x, qmax, tau) for x in t]) / RHOCP
            self.assertLessEqual(np.max(np.abs(err - pred)), dt**2 / 12.0 * (math.pi * qmax / tau) / RHOCP * 1.01)

    def test_meseta_exacta_despues_de_tau(self):
        r, dt, qmax, tau = self.correr(200)
        plateau = qmax * tau / (2 * RHOCP)
        t = r.historial_tiempos_s
        despues = r.historial_puntos_c[t >= tau, 0] - 20.0
        np.testing.assert_allclose(despues, plateau, rtol=1e-11)

    def test_error_decrece_al_refinar_dt(self):
        errores = []
        for n_por_tau in (50, 100, 200):
            r, dt, qmax, tau = self.correr(n_por_tau, horizonte_en_tau=0.5)
            errores.append(np.max(np.abs((r.historial_puntos_c[:, 0] - 20.0) - self.exacto(r.historial_tiempos_s, qmax, tau))))
        self.assertGreater(errores[0] / errores[1], 1.8)     # primer orden
        self.assertGreater(errores[1] / errores[2], 1.8)


class TestScriptDelEjercicio(unittest.TestCase):
    def test_fuente_q_y_linealidad_del_script(self):
        ej = _cargar("ejercicio_fuente_variable_2d", "ejercicio_fuente_variable_2d.py")
        tau = ej.DURACION_FUENTE_S
        self.assertEqual(ej.fuente_q(0.0, 10.0), 0.0)
        self.assertAlmostEqual(ej.fuente_q(tau / 2, 10.0), 10.0)
        self.assertAlmostEqual(ej.fuente_q(tau, 10.0), 0.0, places=12)
        self.assertEqual(ej.fuente_q(tau * 1.0000001, 10.0), 0.0)
        self.assertEqual(ej.fuente_q(-1.0, 10.0), 0.0)
        self.assertAlmostEqual(ej.qmax_por_linealidad(100.0, 25.0, 70.0, 20.0), 200.0)
        with self.assertRaises(ValueError):
            ej.qmax_por_linealidad(100.0, 0.0)


if __name__ == "__main__":
    unittest.main()
