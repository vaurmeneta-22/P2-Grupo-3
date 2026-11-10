"""Pruebas de la secuencia de construcción y de la coherencia de las curvas con la simulación."""
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.analisis import eventos_activacion
from mcplan.thermal.construccion import conjuntos_acumulados, nombre_conjunto
from mcplan.thermal.diffusion_nd import robin
from mcplan.thermal.extremos import SeguimientoExtremos
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.vertidas import simular_vertidas

H = 3600.0


def vertida(identificador, inicio_h, temperatura_c):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    triangulos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", "gid", inicio_h,
                      temperatura_c, "M1", vertices, triangulos, {})


def grilla_tres_capas():
    ids = np.zeros((3, 3, 3), dtype=np.int32)
    for k in range(3):
        ids[:, :, k] = k + 1
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 9.0, 2: 9.0, 3: 9.0})


def simular(grilla, vertidas, seguimiento=None, horas=30.0):
    return simular_vertidas(
        grilla=grilla, vertidas=vertidas, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
        calor_especifico_j_kg_k=900.0, hidratacion=ModeloPulsoCosenoidal(500.0, 12 * H),
        tiempo_final_s=horas * H, paso_tiempo_max_s=1800.0, tiempos_muestreo_s=[4 * H, 20 * H],
        fronteras={cara: robin(8.0, 20.0) for cara in ("x-", "x+", "y-", "y+", "z-", "z+")},
        seguimiento=seguimiento)


class SecuenciaTest(unittest.TestCase):
    def test_el_orden_sigue_las_horas_y_no_los_ids(self):
        # La vertida 3 se coloca primero y la 1 al final.
        grilla = grilla_tres_capas()
        vertidas = [vertida(1, 12.0, 18.0), vertida(2, 6.0, 25.0), vertida(3, 0.0, 30.0)]
        eventos = eventos_activacion(simular(grilla, vertidas, horas=20.0))
        self.assertEqual(list(eventos), [(0.0, (3,)), (6.0, (2,)), (12.0, (1,))])
        self.assertEqual(conjuntos_acumulados(eventos), [(3,), (3, 2), (3, 2, 1)])
        self.assertEqual(nombre_conjunto((3, 2, 1)), "V3+V2+V1")

    def test_simultaneas_se_agrupan_y_se_omiten_las_posteriores_al_horizonte(self):
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0),
                             np.array([[[1, 2, 3, 4]]], dtype=np.int32), {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0})
        vertidas = [vertida(1, 0.0, 18.0), vertida(2, 6.0, 25.0), vertida(3, 6.0, 25.0),
                    vertida(4, 100.0, 25.0)]
        eventos = eventos_activacion(simular(grilla, vertidas, horas=24.0))
        self.assertEqual(list(eventos), [(0.0, (1,)), (6.0, (2, 3))])


class CoherenciaConSimulacionTest(unittest.TestCase):
    def setUp(self):
        self.grilla = grilla_tres_capas()
        # Orden de colocación distinto del orden de IDs: 2 → 3 → 1.
        self.vertidas = [vertida(1, 16.0, 15.0), vertida(2, 0.0, 30.0), vertida(3, 8.0, 18.0)]
        self.extremos = SeguimientoExtremos(self.grilla)
        self.resultado = simular(self.grilla, self.vertidas, self.extremos)

    def fila_en(self, tiempo_s):
        n = int(np.argmin(np.abs(np.asarray(self.extremos.tiempos_s) - tiempo_s)))
        self.assertAlmostEqual(self.extremos.tiempos_s[n], tiempo_s, places=6)
        return n

    def test_curvas_y_campos_guardados_usan_el_mismo_dominio_activo(self):
        s = self.extremos.serie("global")
        for i, tiempo in enumerate(self.resultado.tiempos_s):
            n = self.fila_en(tiempo)
            activas = self.resultado.mascaras_activas[i]
            celdas = self.resultado.temperaturas_c[i][activas]
            self.assertEqual(int(s["n_celdas"][n]), int(activas.sum()))
            self.assertAlmostEqual(s["tmin"][n], celdas.min(), places=9)
            self.assertAlmostEqual(s["tmax"][n], celdas.max(), places=9)
            self.assertAlmostEqual(s["tmedia"][n], celdas.mean(), places=9)

    def test_media_por_vertida_y_exclusion_de_hormigon_no_colocado(self):
        for identificador in (1, 2, 3):
            s = self.extremos.serie(f"vertida_{identificador}")
            horas = s["tiempos_s"] / H
            inicio = {1: 16.0, 2: 0.0, 3: 8.0}[identificador]
            self.assertTrue(np.all(np.isnan(s["tmedia"][horas < inicio])))   # ausente, no 0 °C
            self.assertTrue(np.all(np.isfinite(s["tmedia"][horas >= inicio])))
        # El conjunto activo crece en el orden configurado: 9 → 18 → 27 celdas.
        s = self.extremos.serie("global")
        horas = s["tiempos_s"] / H
        self.assertTrue(np.all(s["n_celdas"][horas < 8.0] == 9))
        self.assertTrue(np.all(s["n_celdas"][(horas >= 8.0) & (horas < 16.0)] == 18))
        self.assertTrue(np.all(s["n_celdas"][horas >= 16.0] == 27))
        self.assertTrue(np.all(s["tmin"] > 0.0))                             # nunca aparece un 0 °C inactivo

    def test_entrada_de_vertida_fria_cambia_el_conjunto_sin_suavizar(self):
        s = self.extremos.serie("global")
        horas = s["tiempos_s"] / H
        n = int(np.flatnonzero(np.isclose(horas, 16.0))[0])
        self.assertEqual(s["tmin"][n], 15.0)                                 # la vertida fría entra a su Tcol
        self.assertLess(s["tmedia"][n], s["tmedia"][n - 1])                  # la media cae al entrar
        self.assertEqual(np.count_nonzero(np.isclose(horas, 16.0)), 1)


if __name__ == "__main__":
    unittest.main()
