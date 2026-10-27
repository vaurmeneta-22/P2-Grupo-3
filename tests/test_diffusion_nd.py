"""Verificación física y numérica del solver explícito 2D/3D."""
import math
import unittest

import numpy as np

from mcplan.thermal.diffusion_nd import (
    DifusionCartesiana, aislado, dirichlet, interpolar_en_punto, neumann, robin)
from mcplan.thermal.verificacion import (
    solucion_sinusoidal, estimar_error_discretizacion, verificar_fuente_uniforme,
    verificar_sinusoidal)

# ---- Parámetros modificables de los casos E05 (unidades en el comentario) ----
E05_T0_C, E05_AMPLITUD_C = 20.0, 30.0                 # °C
E05_2D = dict(forma=(30, 20), longitudes_m=(0.6, 0.4), tiempo_final_s=6 * 3600.0)   # celdas, m, s
E05_3D = dict(forma=(20, 16, 12), longitudes_m=(0.5, 0.4, 0.3), tiempo_final_s=3 * 3600.0)
E05_FUENTE = dict(forma=(6, 5, 4), dx_m=(0.1, 0.1, 0.1), T_inicial_c=20.0,
                  q0_w_m3=500.0, tiempo_final_s=6 * 3600.0)   # celdas, m, °C, W/m³, s

K, RHO, CP = 2.0, 2400.0, 900.0
ALFA = K / (RHO * CP)


def sistema(forma, dx, fronteras=None, k=K, fuente=0.0):
    return DifusionCartesiana(forma=forma, dx_m=dx, conductividad_w_m_k=k,
                              densidad_kg_m3=RHO, calor_especifico_j_kg_k=CP,
                              fuente_w_m3=fuente, fronteras=fronteras)


def todas(condicion, ndim):
    return {c: condicion for c in (("x-", "x+", "y-", "y+", "z-", "z+")[:2 * ndim])}


def modo_seno(forma, L, T0=0.0, A=1.0, t=0.0):
    """Solución exacta (T0 + A·sin·sin[·sin]·exp(−αλt)) en centros; L escalar o tupla por eje."""
    longitudes = (L,) * len(forma) if np.isscalar(L) else tuple(L)
    return solucion_sinusoidal(forma, longitudes, T0, A, ALFA, t)


class TestEquilibrioYConservacion(unittest.TestCase):
    def test_equilibrio_uniforme_2d_3d(self):
        for forma, dx in (((6, 5), (0.1, 0.1)), ((4, 5, 6), (0.1, 0.1, 0.1))):
            s = sistema(forma, dx, todas(aislado(), len(forma)))
            r = s.simular(temperatura_inicial_c=25.0, tiempo_final_s=3600.0,
                          paso_tiempo_max_s=s.paso_maximo_estable_s())
            np.testing.assert_allclose(r.temperatura_c, 25.0, atol=1e-12)

    def test_equilibrio_uniforme_con_dirichlet_igual(self):
        s = sistema((5, 5, 5), (0.1,) * 3, todas(dirichlet(25.0), 3))
        r = s.simular(temperatura_inicial_c=25.0, tiempo_final_s=600.0,
                      paso_tiempo_max_s=s.paso_maximo_estable_s())
        np.testing.assert_allclose(r.temperatura_c, 25.0, atol=1e-12)

    def test_conservacion_energia_adiabatica(self):
        rng = np.random.default_rng(1)
        for forma, dx in (((8, 7), (0.1, 0.12)), ((5, 6, 4), (0.1, 0.1, 0.15))):
            s = sistema(forma, dx, todas(aislado(), len(forma)))
            T0 = rng.uniform(0, 50, forma)
            r = s.simular(temperatura_inicial_c=T0,
                          tiempo_final_s=5 * 3600.0,
                          paso_tiempo_max_s=s.paso_maximo_estable_s())
            self.assertAlmostEqual(r.energia_final_j, r.energia_inicial_j,
                                   delta=1e-9 * abs(r.energia_inicial_j))
            # y el campo se suaviza hacia la media
            self.assertLess(np.std(r.temperatura_c), 0.5 * np.std(T0))

    def test_fuente_uniforme_aislada(self):
        Q, t = 500.0, 6 * 3600.0     # W/m³, s
        for forma, dx in (((5, 5), (0.1, 0.1)), ((4, 4, 4), (0.1,) * 3)):
            s = sistema(forma, dx, todas(aislado(), len(forma)), fuente=Q)
            r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=t,
                          paso_tiempo_max_s=s.paso_maximo_estable_s())
            esperado = 20.0 + Q * t / (RHO * CP)          # T0 + Qdot t/(rho cp)
            np.testing.assert_allclose(r.temperatura_c, esperado, rtol=1e-12)
            self.assertAlmostEqual(r.error_balance_j, 0.0,
                                   delta=1e-9 * r.energia_final_j)

    def test_fuente_uniforme_aislada_temperatura_y_energia(self):
        """T = T_ini + q0 t/(rho cp)  y  E_final − E_inicial = q0·V_total·t."""
        res = verificar_fuente_uniforme(k_w_m_k=K, rho_kg_m3=RHO, cp_j_kg_k=CP, **E05_FUENTE)
        self.assertTrue(res.aprobado, res)
        self.assertLess(res.detalle["error_energia_relativo"], res.detalle["tol_energia_relativa"])

    def test_fuente_cero_equilibrio(self):
        """q0 = 0: T y E no cambian; el criterio de energía no divide entre ΔE_esp = 0."""
        datos = dict(E05_FUENTE, q0_w_m3=0.0)
        res = verificar_fuente_uniforme(k_w_m_k=K, rho_kg_m3=RHO, cp_j_kg_k=CP, **datos)
        self.assertTrue(res.aprobado, res)
        self.assertEqual(res.detalle["dE_esperado_j"], 0.0)
        self.assertLessEqual(res.detalle["error_energia_relativo"], res.detalle["tol_energia_relativa"])
        self.assertLessEqual(abs(res.detalle["dE_j"]), 1e-10 * res.detalle["denominador_energia_j"])

    def test_fuente_por_celda(self):
        Q = np.zeros((4, 4)); Q[1, 1] = 1000.0
        s = sistema((4, 4), (0.1, 0.1), todas(aislado(), 2), fuente=Q)
        r = s.simular(temperatura_inicial_c=0.0, tiempo_final_s=100.0,
                      paso_tiempo_max_s=s.paso_maximo_estable_s())
        self.assertAlmostEqual(r.energia_final_j, 1000.0 * s.volumen_m3 * 100.0, places=6)


class TestSolucionManufacturada(unittest.TestCase):
    """Modo sinusoidal con Dirichlet 0 °C: T = exp(-alfa n pi² t/L²) sin·sin[·sin]."""

    def error(self, forma, L, dt_max, t=3000.0, T0=0.0, A=1.0):
        """Error máx. [°C] vs exacta; L escalar o longitudes por eje; T0, A en °C."""
        ndim = len(forma)
        longitudes = (L,) * ndim if np.isscalar(L) else tuple(L)
        dx = tuple(l / n for l, n in zip(longitudes, forma))
        s = sistema(forma, dx, todas(dirichlet(T0), ndim))
        r = s.simular(temperatura_inicial_c=modo_seno(forma, longitudes, T0, A),
                      tiempo_final_s=t, paso_tiempo_max_s=dt_max(s))
        exacta = modo_seno(forma, longitudes, T0, A, t)
        return float(np.max(np.abs(r.temperatura_c - exacta))), r

    def test_2d_y_3d_cercano_a_la_analitica(self):
        for forma in ((24, 24), (12, 12, 12)):
            e, r = self.error(forma, 0.5, lambda s: 0.9 * s.paso_maximo_estable_s())
            self.assertLess(e, 0.01)     # amplitud inicial 1
            self.assertLess(r.paso_tiempo_s, r.paso_maximo_estable_s)

    def test_refinamiento_conjunto_dx_dt(self):
        """Refinamiento CONJUNTO: dx se divide por 2 y dt por 4 (Fourier fijo).
        El error temporal y el espacial son ambos O(dx²): razón ≈ 4."""
        for ndim, ns in ((2, (16, 32, 64)), (3, (8, 16))):
            errores = []
            for n in ns:
                e, _ = self.error((n,) * ndim, 0.5,
                                  lambda s: 0.4 * s.paso_estable_interior_s())
                errores.append(e)
            for a, b in zip(errores, errores[1:]):
                self.assertGreater(a / b, 3.0, f"{ndim}D errores {errores}")

    def test_convergencia_espacial_sola(self):
        """dt muy pequeño y FIJO: el error restante es espacial (orden 2)."""
        dt = 0.05 * sistema((32, 32), (0.5 / 32,) * 2).paso_maximo_estable_s()
        e = [self.error((n, n), 0.5, lambda s: dt, t=1500.0)[0] for n in (8, 16, 32)]
        self.assertGreater(e[0] / e[1], 3.5)
        self.assertGreater(e[1] / e[2], 3.5)

    def test_convergencia_temporal_sola(self):
        """dx fijo; error contra el mismo esquema espacial con dt mínimo.
        Euler adelantado: orden 1 (razón ≈ 2 al dividir dt por 2)."""
        forma, L, t = (12, 12), 0.5, 3000.0
        s = sistema(forma, (L / 12,) * 2, todas(dirichlet(0.0), 2))
        base = s.paso_maximo_estable_s()
        def corre(dt):
            return s.simular(temperatura_inicial_c=modo_seno(forma, L),
                             tiempo_final_s=t, paso_tiempo_max_s=dt).temperatura_c
        ref = corre(base / 32)
        e = [np.max(np.abs(corre(base / m) - ref)) for m in (1, 2, 4)]
        self.assertGreater(e[0] / e[1], 1.7)
        self.assertGreater(e[1] / e[2], 1.7)


class TestSolucionExactaE05(unittest.TestCase):
    """T0 ≠ 0, A ≠ 1 y lados distintos; Dirichlet T0 en todas las caras, Qdot = 0."""

    error = TestSolucionManufacturada.error      # mismo auxiliar generalizado (T0, A, lados)

    def comprobar(self, caso):
        res = verificar_sinusoidal(T0_c=E05_T0_C, amplitud_c=E05_AMPLITUD_C, k_w_m_k=K,
                                   rho_kg_m3=RHO, cp_j_kg_k=CP, **caso)
        print(f"\n  {res.nombre}: max abs error = {res.error_max_c:.4g} C "
              f"({res.detalle['error_normalizado']:.3g}*|A|), tol = {res.tolerancia_c:.4g} C, "
              f"decay exp(-alfa*lambda*t) = {res.detalle['factor_decaimiento']:.3f}")
        self.assertTrue(res.aprobado, res)
        self.assertLess(res.detalle["factor_decaimiento"], 0.5)     # difusión apreciable
        self.assertLess(res.detalle["dt_s"], res.detalle["dt_limite_s"])
        self.assertLess(abs(res.detalle["T_max_numerica_c"] - res.detalle["T_max_exacta_c"]),
                        res.tolerancia_c)
        return res

    def test_2d(self):
        self.comprobar(E05_2D)

    def test_3d(self):
        self.comprobar(E05_3D)

    def test_varios_pasos_estables_y_dimensiones(self):
        """Distintas dimensiones/lados y fracciones del dt estable. Se exige error <= tol
        (estimación espacial + temporal con el dt efectivo). NO se exige que el error
        baje al reducir solo dt: el error espacial y el de Euler tienen signo opuesto
        y pueden compensarse (cerca de fraccion ~0.6 el error total es mínimo)."""
        casos = (dict(forma=(24, 24), longitudes_m=(0.3, 0.3), tiempo_final_s=3600.0),
                 dict(forma=(36, 12), longitudes_m=(0.9, 0.3), tiempo_final_s=43200.0),
                 dict(forma=(12, 40), longitudes_m=(0.2, 0.8), tiempo_final_s=7200.0),
                 dict(forma=(10, 10, 10), longitudes_m=(0.2, 0.2, 0.2), tiempo_final_s=1800.0),
                 dict(forma=(24, 8, 8), longitudes_m=(0.6, 0.2, 0.2), tiempo_final_s=3600.0),
                 E05_2D, E05_3D)
        for caso in casos:
            for fraccion in (0.2, 0.6, 1.0):
                res = verificar_sinusoidal(T0_c=E05_T0_C, amplitud_c=E05_AMPLITUD_C, k_w_m_k=K,
                                           rho_kg_m3=RHO, cp_j_kg_k=CP, fraccion_paso=fraccion,
                                           **caso)
                self.assertTrue(res.aprobado, (caso, fraccion, res))
                self.assertLessEqual(res.detalle["dt_s"], res.detalle["dt_limite_s"] * (1 + 1e-12))

    def test_estimacion_incluye_error_temporal_con_dt_efectivo(self):
        """A igual malla, un dt mayor aumenta la estimación temporal y no la espacial."""
        args = ((0.6, 0.4), (0.02, 0.02), 30.0, ALFA, 6 * 3600.0)
        e1, t1 = estimar_error_discretizacion(*args, 30.0)
        e2, t2 = estimar_error_discretizacion(*args, 60.0)
        self.assertEqual(e1, e2)
        self.assertAlmostEqual(t2 / t1, 2.0, places=12)

    def test_amplitud_cero_equilibrio(self):
        for caso in (E05_2D, E05_3D):
            res = verificar_sinusoidal(T0_c=E05_T0_C, amplitud_c=0.0, k_w_m_k=K,
                                       rho_kg_m3=RHO, cp_j_kg_k=CP, **caso)
            self.assertTrue(res.aprobado, res)
            self.assertIsNone(res.detalle["error_normalizado"])     # sin dividir entre |A| = 0
            self.assertLessEqual(res.error_max_c, 1e-12 * E05_T0_C)

    def test_error_cae_con_orden_2_al_refinar_conjuntamente(self):
        for ndim, base in ((2, E05_2D), (3, E05_3D)):
            errores = []
            for f in (1, 2):
                forma = tuple(f * n for n in base["forma"])
                res = verificar_sinusoidal(T0_c=E05_T0_C, amplitud_c=E05_AMPLITUD_C, k_w_m_k=K,
                                           rho_kg_m3=RHO, cp_j_kg_k=CP, forma=forma,
                                           longitudes_m=base["longitudes_m"],
                                           tiempo_final_s=base["tiempo_final_s"])
                errores.append(res.error_max_c)
            self.assertGreater(errores[0] / errores[1], 3.5, f"{ndim}D {errores}")

    def test_helper_generalizado_acepta_t0_a_y_lados(self):
        e, _ = self.error((24, 16), (0.6, 0.4), lambda s: 0.9 * s.paso_maximo_estable_s(),
                          t=6 * 3600.0, T0=E05_T0_C, A=E05_AMPLITUD_C)
        self.assertLess(e / E05_AMPLITUD_C, 0.01)


class TestBordes(unittest.TestCase):
    def test_robin_direccion_y_h_cero(self):
        forma, dx = (6, 6), (0.1, 0.1)
        for T0, T_amb, signo in ((60.0, 20.0, -1), (0.0, 20.0, +1)):
            s = sistema(forma, dx, {**todas(aislado(), 2), "x+": robin(10.0, T_amb)})
            r = s.simular(temperatura_inicial_c=T0, tiempo_final_s=3600.0,
                          paso_tiempo_max_s=s.paso_maximo_estable_s())
            self.assertEqual(np.sign(r.temperatura_c.mean() - T0), signo)
            self.assertEqual(np.sign(r.energia_frontera_j["x+"]), signo)
            self.assertAlmostEqual(r.error_balance_j, 0.0, delta=1e-9 * abs(r.energia_final_j))
        s0 = sistema(forma, dx, todas(robin(0.0, 500.0), 2))
        r0 = s0.simular(temperatura_inicial_c=30.0, tiempo_final_s=3600.0,
                        paso_tiempo_max_s=s0.paso_maximo_estable_s())
        np.testing.assert_allclose(r0.temperatura_c, 30.0, atol=1e-12)
        self.assertEqual(r0.energia_frontera_total_j, 0.0)

    def test_robin_estacionario_resistencias_en_serie(self):
        """Dirichlet izq. + Robin der.: q = (T1-Tamb)/(L/k + 1/h)."""
        n, dx, h, T1, Ta = 8, 0.05, 8.0, 50.0, 10.0
        s = sistema((n, 3), (dx, dx), {"x-": dirichlet(T1), "x+": robin(h, Ta),
                                       "y-": aislado(), "y+": aislado()})
        r = s.simular(temperatura_inicial_c=Ta, tiempo_final_s=1.0e6,
                      paso_tiempo_max_s=s.paso_maximo_estable_s())
        q_esperado = (T1 - Ta) / (n * dx / K + 1.0 / h)          # W/m²
        area = s.area_cara_m2[0] * 3                            # 3 celdas en y
        t = r.tiempo_s
        # potencia que sale por Robin al final: G (T - Ta) por celda de borde
        T_borde = r.temperatura_c[-1, :]
        q_robin = float(np.sum(s.G_borde["x+"] * (T_borde - Ta))) / area
        self.assertAlmostEqual(q_robin, q_esperado, delta=1e-4 * q_esperado)

    def test_neumann_no_nulo_signo_y_balance(self):
        q, t = 200.0, 3600.0
        forma, dx = (5, 4, 3), (0.1, 0.1, 0.1)
        for signo in (+1, -1):
            s = sistema(forma, dx, {**todas(aislado(), 3), "z-": neumann(signo * q)})
            r = s.simular(temperatura_inicial_c=20.0, tiempo_final_s=t,
                          paso_tiempo_max_s=s.paso_maximo_estable_s())
            area = dx[0] * dx[1] * forma[0] * forma[1]
            self.assertAlmostEqual(r.energia_frontera_j["z-"], signo * q * area * t, places=6)
            self.assertEqual(np.sign(r.temperatura_c.mean() - 20.0), signo)
            self.assertAlmostEqual(r.energia_final_j - r.energia_inicial_j,
                                   signo * q * area * t, delta=1e-9 * r.energia_final_j)

    def test_neumann_equivale_a_gradiente_1d(self):
        """q_in = -k dT/dx en x mínimo (convenio del solver 1D con gradiente)."""
        s = sistema((4, 3), (0.1, 0.1), {**todas(aislado(), 2), "x-": neumann(-K * (-100.0))})
        # dT/dx = -100 K/m en x=0  ->  calor entra hacia +x: q_in = +200 W/m²
        r = s.simular(temperatura_inicial_c=0.0, tiempo_final_s=10.0,
                      paso_tiempo_max_s=s.paso_maximo_estable_s())
        self.assertGreater(r.energia_frontera_j["x-"], 0.0)


class TestMaterialesYEstabilidad(unittest.TestCase):
    def test_dos_materiales_flujo_continuo_y_resistencia(self):
        n, dx = 10, 0.05
        k = np.empty((n, 2)); k[:5] = 1.0; k[5:] = 4.0          # W/(m K)
        s = sistema((n, 2), (dx, dx), {"x-": dirichlet(100.0), "x+": dirichlet(0.0),
                                       "y-": aislado(), "y+": aislado()}, k=k)
        r = s.simular(temperatura_inicial_c=50.0, tiempo_final_s=1.0e6,
                      paso_tiempo_max_s=s.paso_maximo_estable_s())
        R = 0.25 / 1.0 + 0.25 / 4.0                               # Σ L/k [m² K/W]
        q_esperado = 100.0 / R
        # flujo en cada cara interna (por área): constante e igual al esperado
        f = s.flujos_internos_w(r.temperatura_c, 0)               # (n-1, 2) [W]
        q = -f / s.area_cara_m2[0]
        np.testing.assert_allclose(q, q_esperado, rtol=1e-4)
        # conductividad efectiva de la cara de la interfaz = media armónica
        k_f = 2 * 1.0 * 4.0 / (1.0 + 4.0)
        self.assertAlmostEqual(s.G_interna[0][4, 0] * dx / s.area_cara_m2[0], k_f)

    def test_flujos_iguales_y_opuestos(self):
        rng = np.random.default_rng(3)
        k = rng.uniform(0.5, 3.0, (5, 5, 5))
        s = sistema((5, 5, 5), (0.1,) * 3, todas(aislado(), 3), k=k)
        T = rng.uniform(0, 40, (5, 5, 5))
        p, _ = s._derivada_energia_w(T, 0.0)
        self.assertAlmostEqual(float(p.sum()), 0.0, delta=1e-9)   # fronteras aisladas

    def test_limite_homogeneo_recupera_formula_clasica(self):
        for forma, dx in (((9, 9), (0.1, 0.2)), ((6, 6, 6), (0.1, 0.2, 0.3))):
            s = sistema(forma, dx, todas(aislado(), len(forma)))
            esperado = 0.5 / (ALFA * sum(1 / d**2 for d in dx))
            self.assertAlmostEqual(s.paso_maximo_estable_s(), esperado, delta=1e-9 * esperado)

    def test_frontera_e_interfaz_endurecen_el_limite(self):
        dx = (0.1, 0.1)
        aisl = sistema((6, 6), dx, todas(aislado(), 2))
        dir_ = sistema((6, 6), dx, todas(dirichlet(0.0), 2))
        self.assertLess(dir_.paso_maximo_estable_s(), aisl.paso_maximo_estable_s())
        rob = sistema((6, 6), dx, todas(robin(1e6, 0.0), 2))     # h enorme ~ Dirichlet
        self.assertAlmostEqual(rob.paso_maximo_estable_s(), dir_.paso_maximo_estable_s(), delta=1e-3 * dir_.paso_maximo_estable_s())
        k = np.ones((6, 6)) * 1.0; k[3:] = 50.0                  # interfaz muy contrastada
        mat = sistema((6, 6), dx, todas(aislado(), 2), k=k)
        self.assertLess(mat.paso_maximo_estable_s(), aisl.paso_maximo_estable_s())

    def test_rechaza_paso_inestable_y_acepta_el_limite(self):
        s = sistema((6, 6, 6), (0.1,) * 3, todas(dirichlet(0.0), 3))
        lim = s.paso_maximo_estable_s()
        with self.assertRaisesRegex(ValueError, "inestable"):
            s.simular(temperatura_inicial_c=1.0, tiempo_final_s=3 * lim, paso_tiempo_max_s=1.5 * lim)
        r = s.simular(temperatura_inicial_c=1.0, tiempo_final_s=100 * lim, paso_tiempo_max_s=lim)
        self.assertTrue(np.all(np.isfinite(r.temperatura_c)))
        self.assertLessEqual(np.max(np.abs(r.temperatura_c)), 1.0 + 1e-12)   # principio del máximo

    def test_entradas_invalidas(self):
        with self.assertRaises(ValueError):
            sistema((4, 4), (0.1, 0.1), k=0.0)
        with self.assertRaises(ValueError):
            sistema((4, 4), (0.1, -0.1))
        with self.assertRaises(ValueError):
            sistema((4, 4), (0.1, 0.1), {"z+": aislado()})
        with self.assertRaises(ValueError):
            robin(-1.0, 20.0)


class TestInterpolacion(unittest.TestCase):
    def test_punto_fisico_igual_en_mallas_distintas(self):
        for ndim in (2, 3):
            valores = []
            for n in (6, 12):
                dx = tuple(1.0 / n for _ in range(ndim))
                ejes = np.meshgrid(*[(np.arange(n) + .5) / n for _ in range(ndim)], indexing="ij")
                T = 3.0 + 2.0 * ejes[0] - ejes[1]            # campo lineal: exacto
                valores.append(interpolar_en_punto(T, dx, (0.5,) * ndim))
            self.assertAlmostEqual(valores[0], 3.0 + 1.0 - 0.5, places=12)
            self.assertAlmostEqual(valores[0], valores[1], places=12)


if __name__ == "__main__":
    unittest.main()
