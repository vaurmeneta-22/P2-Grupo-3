"""Pruebas de difusión sobre celdas activas de una grilla 3D irregular."""
import unittest

import numpy as np

from mcplan.thermal.diffusion_nd import DifusionCartesiana, aislado, robin


def crear_solver(mascara, fronteras=None):
    return DifusionCartesiana(
        forma=mascara.shape,
        dx_m=(1.0, 1.0, 1.0),
        conductividad_w_m_k=1.0,
        densidad_kg_m3=1.0,
        calor_especifico_j_kg_k=1.0,
        fronteras=fronteras,
        mascara_activa=mascara,
    )


class DifusionConMascaraTests(unittest.TestCase):
    def test_dos_voxeles_activos_intercambian_calor_y_conservan_energia(self):
        mascara = np.ones((2, 1, 1), dtype=bool)
        solver = crear_solver(mascara)
        inicial = np.array([0.0, 100.0]).reshape(2, 1, 1)
        resultado = solver.simular(temperatura_inicial_c=inicial,
                                   tiempo_final_s=2.0, paso_tiempo_max_s=0.5)
        self.assertGreater(resultado.temperatura_c[0, 0, 0], inicial[0, 0, 0])
        self.assertLess(resultado.temperatura_c[1, 0, 0], inicial[1, 0, 0])
        self.assertAlmostEqual(resultado.error_balance_j, 0.0, places=12)

    def test_vacio_separa_voxeles_y_no_acumula_energia(self):
        mascara = np.array([True, False, True]).reshape(3, 1, 1)
        solver = crear_solver(mascara)
        inicial = np.array([0.0, 999.0, 100.0]).reshape(3, 1, 1)
        resultado = solver.simular(temperatura_inicial_c=inicial,
                                   tiempo_final_s=1.0, paso_tiempo_max_s=0.2)
        np.testing.assert_allclose(resultado.temperatura_c.ravel(), [0.0, 999.0, 100.0])
        self.assertAlmostEqual(resultado.energia_inicial_j, 100.0)
        self.assertAlmostEqual(resultado.energia_final_j, 100.0)

    def test_robin_actua_tambien_en_superficies_hacia_voxel_vacio(self):
        mascara = np.array([False, True, False]).reshape(3, 1, 1)
        fronteras = {"x-": robin(1.0, 0.0), "x+": robin(1.0, 0.0),
                     "y-": aislado(), "y+": aislado(),
                     "z-": aislado(), "z+": aislado()}
        solver = crear_solver(mascara, fronteras)
        inicial = np.array([0.0, 10.0, 0.0]).reshape(3, 1, 1)
        resultado = solver.simular(temperatura_inicial_c=inicial,
                                   tiempo_final_s=0.1, paso_tiempo_max_s=0.1)
        self.assertLess(resultado.temperatura_c[1, 0, 0], 10.0)
        self.assertLess(resultado.energia_frontera_j["x-"], 0.0)
        self.assertLess(resultado.energia_frontera_j["x+"], 0.0)
        self.assertAlmostEqual(resultado.error_balance_j, 0.0, places=12)

    def test_mascara_puede_actualizarse_al_activar_otra_vertida(self):
        inicial = np.array([True, False]).reshape(2, 1, 1)
        solver = crear_solver(inicial)
        solver.actualizar_mascara_activa(np.ones((2, 1, 1), dtype=bool))
        self.assertEqual(np.count_nonzero(solver.mascara_activa), 2)
        self.assertGreater(np.count_nonzero(solver.G_interna[0]), 0)

    def test_rechaza_mascara_con_forma_incorrecta(self):
        with self.assertRaises(ValueError):
            crear_solver(np.ones((2, 1, 1), dtype=bool)).actualizar_mascara_activa(
                np.ones((1, 1, 1), dtype=bool))


if __name__ == "__main__":
    unittest.main(verbosity=2)
