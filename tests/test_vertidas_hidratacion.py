"""Pruebas de fuente adiabática y activación temporal por vertida."""
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.hidratacion import (
    ModeloHidratacionExponencial,
    ModeloPulsoCosenoidal,
)
from mcplan.thermal.vertidas import simular_vertidas


def crear_vertida(identificador, inicio_h, temperatura_c):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                         [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    caras = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                      [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                      [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], dtype=int)
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", f"guid-{identificador}",
                      inicio_h, temperatura_c, "M1", vertices, caras,
                      {"VertidaID": identificador, "Inicio_h": inicio_h,
                       "TemperaturaColocacion_C": temperatura_c})


class HidratacionTests(unittest.TestCase):
    def test_fuente_es_cero_antes_de_colocar_y_deriva_de_curva_adiabatica(self):
        modelo = ModeloHidratacionExponencial(20.0, 3600.0)
        rho, cp = 2400.0, 900.0
        q = modelo.generacion_w_m3(np.array([-1.0, 0.0]), rho, cp)
        self.assertEqual(q[0], 0.0)
        self.assertAlmostEqual(q[1], rho * cp * 20.0 / 3600.0)
        self.assertAlmostEqual(modelo.temperatura_adiabatica_c(18.0, 0.0), 18.0)

    def test_incremento_adiabatico_converge_al_delta_configurado(self):
        modelo = ModeloHidratacionExponencial(18.0, 7200.0)
        rho, cp = 2200.0, 850.0
        tiempos = np.linspace(0.0, 10 * modelo.tiempo_caracteristico_s, 10001)
        q = modelo.generacion_w_m3(tiempos, rho, cp)
        energia = float(np.trapezoid(q, tiempos))
        self.assertAlmostEqual(energia / (rho * cp), 18.0, delta=0.01)

    def test_dos_vertidas_se_activan_con_tiempo_y_temperatura_de_sus_inputs(self):
        ids = np.array([1, 2], dtype=np.int32).reshape(2, 1, 1)
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids,
                             {1: 1.0, 2: 1.0})
        vertidas = [crear_vertida(1, 0.0, 18.0), crear_vertida(2, 0.5, 25.0)]
        resultado = simular_vertidas(
            grilla=grilla, vertidas=vertidas,
            conductividad_w_m_k=1.0, densidad_kg_m3=1000.0,
            calor_especifico_j_kg_k=1000.0,
            hidratacion=ModeloHidratacionExponencial(5.0, 3600.0),
            tiempo_final_s=3600.0, paso_tiempo_max_s=60.0,
        )
        self.assertEqual(resultado.ids_activados, ((1,), (2,), ()))
        self.assertEqual([int(m.sum()) for m in resultado.mascaras_activas], [1, 2, 2])
        self.assertAlmostEqual(resultado.temperaturas_c[1][1, 0, 0], 25.0)
        self.assertGreater(resultado.temperaturas_c[-1][1, 0, 0], 25.0)
        self.assertAlmostEqual(resultado.error_balance_j, 0.0, delta=1e-5)


class FuenteCosenoidalTests(unittest.TestCase):
    def setUp(self):
        self.duracion_s = 48.0 * 3600.0
        self.modelo = ModeloPulsoCosenoidal(500.0, self.duracion_s)

    def test_pulso_empieza_en_cero_alcanza_qmax_y_se_apaga(self):
        edades = np.array([
            -1.0,
            0.0,
            self.duracion_s / 4.0,
            self.duracion_s / 2.0,
            self.duracion_s,
            self.duracion_s + 1.0,
        ])
        generacion = self.modelo.generacion_w_m3(edades, 2400.0, 900.0)
        np.testing.assert_allclose(generacion, [0.0, 0.0, 250.0, 500.0, 0.0, 0.0])

    def test_pulso_entrega_20_c_adiabaticos_con_los_parametros_del_proyecto(self):
        temperatura = self.modelo.temperatura_adiabatica_c(
            18.0, self.duracion_s, 2400.0, 900.0
        )
        self.assertAlmostEqual(temperatura, 38.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
