"""Pruebas del registro de hidratación y flujos (E11): coincidencia con el solver y figuras."""
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.bordes import ModeloBordes, ModeloSuelo, Proteccion
from mcplan.thermal.diffusion_nd import robin
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.registro_flujos import RegistroEnergetico, mascara_base
from mcplan.thermal.vertidas import simular_vertidas
from mcplan.visualization import figuras_e11 as fig11

H = 3600.0
CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")


def vertida(identificador, inicio_h, temperatura_c):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    triangulos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", "gid", inicio_h, temperatura_c, "M1",
                      vertices, triangulos, {})


def grilla_dos_capas():
    ids = np.zeros((3, 3, 2), dtype=np.int32)
    ids[:, :, 0] = 1
    ids[:, :, 1] = 2
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 9.0, 2: 9.0})


def correr(observadores=None, *, aire_c=10.0, retiro_h=12.0, temperatura_c=(40.0, 35.0), horas=30.0):
    grilla = grilla_dos_capas()
    bordes = ModeloBordes(
        grilla, robin(10.0, aire_c), suelo=ModeloSuelo(8.0, 1.5, 0.75),
        protecciones=(Proteccion(1, "moldaje", ("x+", "x-", "y+", "y-"), 2.0, retiro_h * H),))
    hidratacion = ModeloPulsoCosenoidal(400.0, 20 * H)
    resultado = simular_vertidas(
        grilla=grilla, vertidas=[vertida(1, 0.0, temperatura_c[0]), vertida(2, 6.0, temperatura_c[1])],
        conductividad_w_m_k=2.0, densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
        hidratacion=hidratacion, tiempo_final_s=horas * H, paso_tiempo_max_s=1800.0,
        tiempos_muestreo_s=[5 * H, 10 * H, 20 * H], bordes=bordes, seguimiento=observadores)
    return grilla, bordes, hidratacion, resultado


class RegistroTest(unittest.TestCase):
    def registro(self, grilla):
        return RegistroEnergetico(grilla, caras=[(1, "x+"), (2, "x+")], celdas_suelo=mascara_base(grilla))

    def test_los_registros_coinciden_con_la_fuente_y_los_flujos_del_solver(self):
        grilla = grilla_dos_capas()
        reg = self.registro(grilla)
        _, _, hidratacion, resultado = correr([reg])
        # Energía generada y energía por la frontera: las mismas que acumula el solver.
        self.assertAlmostEqual(reg.energia_fuente_j / resultado.energia_fuente_j, 1.0, places=12)
        self.assertAlmostEqual(reg.energia_frontera_j / resultado.energia_frontera_j, 1.0, places=10)
        # Qdot registrado = modelo de hidratación aplicado a la edad de cada vertida.
        t = np.asarray(reg.t_s)
        for i, inicio_h in ((1, 0.0), (2, 6.0)):
            esperado = hidratacion.generacion_w_m3(t - inicio_h * H, 2400.0, 900.0)
            registrado = np.asarray(reg.qdot_w_m3[i])
            validos = np.isfinite(registrado)
            self.assertTrue(np.allclose(registrado[validos], esperado[validos], rtol=0, atol=1e-9))
            self.assertTrue(np.all(np.isnan(registrado[t < inicio_h * H - 1e-9])))   # ausente antes de colocarse
        # Potencia total = Σ potencias por vertida (activas).
        suma = np.nansum([reg.potencia_vertida_w[i] for i in reg.ids], axis=0)
        self.assertTrue(np.allclose(reg.potencia_total_w, suma))

    def test_la_solucion_termica_no_cambia_con_los_observadores(self):
        grilla = grilla_dos_capas()
        _, _, _, con = correr([self.registro(grilla)])
        _, _, _, sin = correr(None)
        for a, b in zip(con.temperaturas_c, sin.temperaturas_c):
            self.assertTrue(np.array_equal(a, b))
        self.assertEqual(con.energia_fuente_j, sin.energia_fuente_j)
        self.assertEqual(con.energia_frontera_j, sin.energia_frontera_j)

    def test_convenio_de_signo_y_unidades_w_vs_w_m2(self):
        grilla = grilla_dos_capas()
        caliente = self.registro(grilla)
        correr([caliente], aire_c=5.0, temperatura_c=(40.0, 40.0))
        self.assertLess(caliente.cara_w_m2[(1, "x+")][1], 0.0)        # hormigón caliente: el calor SALE
        self.assertLess(caliente.suelo_w_m2[1], 0.0)
        frio = self.registro(grilla)
        correr([frio], aire_c=35.0, temperatura_c=(5.0, 5.0))
        self.assertGreater(frio.cara_w_m2[(1, "x+")][1], 0.0)          # hormigón frío: el calor ENTRA
        # W/m² = W / (celdas · área); el suelo usa la cara z- (1 m²) de 9 celdas.
        n = caliente.suelo_n[3]
        self.assertEqual(n, 9)
        self.assertAlmostEqual(caliente.suelo_w_m2[3], caliente.suelo_w[3] / (9 * 1.0))
        n_cara = caliente.cara_n[(1, "x+")][3]
        self.assertAlmostEqual(caliente.cara_w_m2[(1, "x+")][3], caliente.cara_w[(1, "x+")][3] / (n_cara * 1.0))

    def test_el_flujo_salta_al_retirar_la_proteccion(self):
        grilla = grilla_dos_capas()
        reg = self.registro(grilla)
        _, _, _, resultado = correr([reg], retiro_h=12.0, aire_c=5.0, temperatura_c=(40.0, 40.0))
        horas = reg.horas
        antes = reg.cara_w_m2[(1, "x+")][int(np.flatnonzero(horas < 12.0)[-1])]
        despues = reg.cara_w_m2[(1, "x+")][int(np.flatnonzero(horas >= 12.0)[0])]
        self.assertLess(despues, antes - 5.0)                           # sale mucho más calor sin moldaje
        r = np.asarray(reg.cara_resistencia[(1, "x+")])
        self.assertTrue(np.allclose(r[horas < 12.0], 2.0))              # protección vigente
        self.assertTrue(np.allclose(r[horas >= 12.0], 0.0))             # retirada
        self.assertEqual(len(resultado.retiros), 1)
        # La vertida 2 no tiene protección: su resistencia es siempre cero.
        r2 = np.asarray(reg.cara_resistencia[(2, "x+")])
        self.assertTrue(np.all(r2[np.isfinite(r2)] == 0.0))

    def test_suelo_en_serie_separada_de_las_caras_laterales(self):
        grilla = grilla_dos_capas()
        reg = self.registro(grilla)
        correr([reg], aire_c=5.0)
        total_lateral_mas_suelo = np.asarray(reg.frontera_total_w)
        self.assertTrue(np.all(np.isfinite(reg.suelo_w)))
        # El total de la frontera incluye el suelo, y no es igual a él (hay caras de aire).
        self.assertFalse(np.allclose(total_lateral_mas_suelo, reg.suelo_w))

    def test_figuras_y_csv(self):
        grilla = grilla_dos_capas()
        reg = self.registro(grilla)
        _, bordes, _, resultado = correr([reg])
        datos = {"activaciones": [(0.0, (1,)), (6.0, (2,))], "retiros": [(12.0, ("V1:moldaje",))],
                 "protecciones": list(bordes.protecciones), "inicio_h": {1: 0.0, 2: 6.0},
                 "aire_horas": np.arange(0, 31.0), "aire_c": [10.0] * 31}
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            f1 = fig11.figura_hidratacion(reg, datos, 30.0, tmp / "h.png")
            f2 = fig11.figura_flujos(reg, datos, 30.0, tmp / "f.png")
            f3 = fig11.figura_por_edad(reg, datos, tmp / "e.png")
            self.assertEqual(f1._suptitle.get_text(), "Generación de calor por hidratación")
            self.assertIn("positivo = entra", f2._suptitle.get_text())
            for nombre in ("h.png", "f.png", "e.png"):
                self.assertGreater((tmp / nombre).stat().st_size, 5000)
            reg.guardar_hidratacion_csv(tmp / "h.csv")
            reg.guardar_flujos_csv(tmp / "f.csv")
            lineas = (tmp / "f.csv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(lineas), len(reg.t_s) + 1)
            self.assertIn("suelo_W_m2", lineas[0])
            self.assertIn("V1_x+_W_m2", lineas[0])
            self.assertEqual(len((tmp / "h.csv").read_text(encoding="utf-8").splitlines()), len(reg.t_s) + 1)
        for figura in (f1, f2, f3):
            plt.close(figura)


if __name__ == "__main__":
    unittest.main()
