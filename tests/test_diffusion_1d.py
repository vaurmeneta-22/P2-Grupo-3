import math
import unittest

import numpy as np

from mcplan.thermal.diffusion_1d import (
    calcular_difusividad_termica,
    simular_difusion_1d,
    simular_dirichlet,
)


class PruebasDifusion1D(unittest.TestCase):
    def test_compara_con_solucion_analitica_y_converge(self):
        difusividad = 1.0e-4
        longitud = 1.0
        temperatura_base = 20.0
        amplitud = 30.0
        tiempo_final = 5000.0
        errores = []

        for numero_nodos in (21, 41):
            posicion = np.linspace(0.0, longitud, numero_nodos)
            temperatura_inicial = temperatura_base + amplitud * np.sin(
                np.pi * posicion / longitud
            )
            resultado = simular_dirichlet(
                temperatura_inicial_c=temperatura_inicial,
                longitud_m=longitud,
                difusividad_m2_s=difusividad,
                temperatura_izquierda_c=temperatura_base,
                temperatura_derecha_c=temperatura_base,
                tiempo_final_s=tiempo_final,
                paso_tiempo_max_s=0.45 * (posicion[1] - posicion[0]) ** 2 / difusividad,
            )
            temperatura_exacta = temperatura_base + amplitud * np.sin(
                np.pi * posicion / longitud
            ) * math.exp(-difusividad * (np.pi / longitud) ** 2 * tiempo_final)
            errores.append(np.max(np.abs(resultado.temperatura_c - temperatura_exacta)))
            self.assertAlmostEqual(resultado.posicion_m[-1], longitud)
            self.assertAlmostEqual(
                resultado.paso_tiempo_s * resultado.numero_pasos, tiempo_final
            )
            self.assertLessEqual(resultado.numero_fourier, 0.5)

        self.assertLess(errores[1], errores[0] / 3.0)

    def test_mantiene_las_temperaturas_dirichlet(self):
        resultado = simular_dirichlet(
            temperatura_inicial_c=[10.0, 25.0, 40.0, 25.0, 10.0],
            longitud_m=2.0,
            difusividad_m2_s=1.0e-5,
            temperatura_izquierda_c=5.0,
            temperatura_derecha_c=15.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=10.0,
        )
        self.assertEqual(resultado.temperatura_c[0], 5.0)
        self.assertEqual(resultado.temperatura_c[-1], 15.0)

    def test_equilibrio_uniforme_no_cambia(self):
        resultado = simular_dirichlet(
            temperatura_inicial_c=np.full(5, 20.0),
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            temperatura_izquierda_c=20.0,
            temperatura_derecha_c=20.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=10.0,
        )
        np.testing.assert_allclose(resultado.temperatura_c, 20.0, atol=1e-12)

    def test_rechaza_un_paso_inestable(self):
        with self.assertRaisesRegex(ValueError, "inestable"):
            simular_dirichlet(
                temperatura_inicial_c=[20.0, 30.0, 20.0],
                longitud_m=1.0,
                difusividad_m2_s=1.0,
                temperatura_izquierda_c=20.0,
                temperatura_derecha_c=20.0,
                tiempo_final_s=4.0,
                paso_tiempo_max_s=4.0,
            )

    def test_calcula_difusividad_desde_propiedades(self):
        difusividad = calcular_difusividad_termica(2.0, 2400.0, 900.0)
        self.assertAlmostEqual(difusividad, 2.0 / (2400.0 * 900.0))

    def test_rechaza_una_temperatura_inicial_con_forma_incorrecta(self):
        with self.assertRaisesRegex(ValueError, "arreglo 1D"):
            simular_dirichlet(
                temperatura_inicial_c=[[20.0, 20.0], [20.0, 20.0]],
                longitud_m=1.0,
                difusividad_m2_s=1.0e-5,
                temperatura_izquierda_c=20.0,
                temperatura_derecha_c=20.0,
                tiempo_final_s=10.0,
                paso_tiempo_max_s=1.0,
            )

    def test_neumann_cero_mantiene_equilibrio_uniforme(self):
        resultado = simular_difusion_1d(
            temperatura_inicial_c=np.full(5, 23.0),
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            condicion_izquierda="neumann",
            valor_izquierdo=0.0,
            condicion_derecha="neumann",
            valor_derecho=0.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=10.0,
        )
        np.testing.assert_allclose(resultado.temperatura_c, 23.0, atol=1e-12)

    def test_neumann_cero_conserva_energia_trapezoidal(self):
        temperatura_inicial = np.array([20.0, 35.0, 50.0, 35.0, 20.0])
        resultado = simular_difusion_1d(
            temperatura_inicial_c=temperatura_inicial,
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            condicion_izquierda="neumann",
            valor_izquierdo=0.0,
            condicion_derecha="neumann",
            valor_derecho=0.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=5.0,
        )
        promedio_inicial = np.trapezoid(temperatura_inicial)
        promedio_final = np.trapezoid(resultado.temperatura_c)
        self.assertAlmostEqual(promedio_final, promedio_inicial, places=10)

    def test_permite_combinar_dirichlet_y_neumann(self):
        resultado = simular_difusion_1d(
            temperatura_inicial_c=[25.0, 30.0, 35.0, 40.0, 45.0],
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            condicion_izquierda="dirichlet",
            valor_izquierdo=20.0,
            condicion_derecha="neumann",
            valor_derecho=0.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=5.0,
        )
        self.assertEqual(resultado.temperatura_c[0], 20.0)
        self.assertNotEqual(resultado.temperatura_c[-1], 45.0)

    def test_robin_en_equilibrio_con_el_ambiente_no_cambia(self):
        resultado = simular_difusion_1d(
            temperatura_inicial_c=np.full(5, 20.0),
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            condicion_izquierda="robin",
            valor_izquierdo=None,
            condicion_derecha="robin",
            valor_derecho=None,
            conductividad_w_m_k=2.0,
            coeficiente_conveccion_izquierdo_w_m2_k=10.0,
            coeficiente_conveccion_derecho_w_m2_k=10.0,
            temperatura_ambiente_izquierda_c=20.0,
            temperatura_ambiente_derecha_c=20.0,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=10.0,
        )
        np.testing.assert_allclose(resultado.temperatura_c, 20.0, atol=1e-12)

    def test_robin_enfria_hacia_el_ambiente(self):
        resultado = simular_difusion_1d(
            temperatura_inicial_c=np.full(21, 50.0),
            longitud_m=1.0,
            difusividad_m2_s=1.0e-4,
            condicion_izquierda="robin",
            valor_izquierdo=None,
            condicion_derecha="robin",
            valor_derecho=None,
            conductividad_w_m_k=2.0,
            coeficiente_conveccion_izquierdo_w_m2_k=10.0,
            coeficiente_conveccion_derecho_w_m2_k=10.0,
            temperatura_ambiente_izquierda_c=20.0,
            temperatura_ambiente_derecha_c=20.0,
            tiempo_final_s=5000.0,
            paso_tiempo_max_s=0.35 * 0.05**2 / 1.0e-4,
        )
        self.assertLess(resultado.temperatura_c[0], 50.0)
        self.assertLess(resultado.temperatura_c[-1], 50.0)
        self.assertGreaterEqual(float(np.min(resultado.temperatura_c)), 20.0)

    def test_robin_con_h_cero_equivale_a_neumann_aislado(self):
        parametros = dict(
            temperatura_inicial_c=[20.0, 30.0, 40.0, 30.0, 20.0],
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            tiempo_final_s=100.0,
            paso_tiempo_max_s=5.0,
        )
        resultado_neumann = simular_difusion_1d(
            **parametros,
            condicion_izquierda="neumann",
            valor_izquierdo=0.0,
            condicion_derecha="neumann",
            valor_derecho=0.0,
        )
        resultado_robin = simular_difusion_1d(
            **parametros,
            condicion_izquierda="robin",
            valor_izquierdo=None,
            condicion_derecha="robin",
            valor_derecho=None,
            conductividad_w_m_k=2.0,
            coeficiente_conveccion_izquierdo_w_m2_k=0.0,
            coeficiente_conveccion_derecho_w_m2_k=0.0,
            temperatura_ambiente_izquierda_c=20.0,
            temperatura_ambiente_derecha_c=20.0,
        )
        np.testing.assert_allclose(
            resultado_robin.temperatura_c, resultado_neumann.temperatura_c
        )

    def test_robin_reduce_el_limite_de_estabilidad(self):
        with self.assertRaisesRegex(ValueError, "inestable"):
            simular_difusion_1d(
                temperatura_inicial_c=[50.0, 50.0, 50.0],
                longitud_m=1.0,
                difusividad_m2_s=1.0e-3,
                condicion_izquierda="robin",
                valor_izquierdo=None,
                condicion_derecha="neumann",
                valor_derecho=0.0,
                conductividad_w_m_k=1.0,
                coeficiente_conveccion_izquierdo_w_m2_k=100.0,
                temperatura_ambiente_izquierda_c=20.0,
                tiempo_final_s=125.0,
                paso_tiempo_max_s=125.0,
            )

    def test_robin_acepta_ambiente_variable_en_el_tiempo(self):
        resultado = simular_difusion_1d(
            temperatura_inicial_c=np.full(5, 20.0),
            longitud_m=1.0,
            difusividad_m2_s=1.0e-5,
            condicion_izquierda="robin",
            valor_izquierdo=None,
            condicion_derecha="neumann",
            valor_derecho=0.0,
            conductividad_w_m_k=2.0,
            coeficiente_conveccion_izquierdo_w_m2_k=10.0,
            temperatura_ambiente_izquierda_c=lambda tiempo_s: (
                20.0 + 10.0 * math.sin(math.pi * tiempo_s / 2000.0)
            ),
            tiempo_final_s=1000.0,
            paso_tiempo_max_s=1.0,
        )
        self.assertGreater(resultado.temperatura_c[0], 20.0)


if __name__ == "__main__":
    unittest.main()
