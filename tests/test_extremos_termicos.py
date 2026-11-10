"""Pruebas de los extremos térmicos: ΔT entre celdas y contraste con caras exteriores."""
import csv
import itertools
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.diffusion_nd import DifusionCartesiana, dirichlet, neumann, robin
from mcplan.thermal.extremos import SeguimientoExtremos
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.seguimiento import CARAS, SeguimientoPuntos, resolver_puntos, temperatura_cara_estimada
from mcplan.thermal.vertidas import simular_vertidas

H = 3600.0


def vertida(identificador, inicio_h, temperatura_c):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    triangulos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", "gid", inicio_h,
                      temperatura_c, "M1", vertices, triangulos, {})


def grilla_dos_capas():
    """Columna 3x3 en planta: vertida 1 en z=0 y vertida 2 encima, dx = 1 m."""
    ids = np.zeros((3, 3, 2), dtype=np.int32)
    ids[:, :, 0], ids[:, :, 1] = 1, 2
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 9.0, 2: 9.0})


def fronteras(condicion):
    return {cara: condicion for cara in CARAS}


def simular(grilla, vertidas, seguimiento=None, horas=24.0, muestras=(8 * H,)):
    return simular_vertidas(
        grilla=grilla, vertidas=vertidas, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
        calor_especifico_j_kg_k=900.0, hidratacion=ModeloPulsoCosenoidal(500.0, 12 * H),
        tiempo_final_s=horas * H, paso_tiempo_max_s=1800.0, tiempos_muestreo_s=list(muestras),
        fronteras=fronteras(robin(8.0, 20.0)), seguimiento=seguimiento)


def solver_linea(mascara, condicion, n=4):
    """Línea de n celdas de 1 m (x), con las 6 caras exteriores bajo ``condicion``."""
    ids = np.zeros((n, 1, 1), dtype=np.int32)
    ids[np.asarray(mascara, dtype=bool), 0, 0] = 1
    grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: float(ids.sum())})
    solver = DifusionCartesiana(
        forma=(n, 1, 1), dx_m=(1.0, 1.0, 1.0), conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
        calor_especifico_j_kg_k=900.0, fronteras=fronteras(condicion),
        mascara_activa=np.asarray(mascara, dtype=bool).reshape(n, 1, 1))
    return grilla, solver


class CampoUniformeTest(unittest.TestCase):
    def test_campo_uniforme_no_tiene_diferencia_entre_celdas(self):
        grilla, solver = solver_linea([1, 1, 1, 1], robin(8.0, 20.0))
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.full((4, 1, 1), 30.0), solver)
        s = extremos.serie("global")
        self.assertEqual(s["tmin"][0], 30.0)
        self.assertEqual(s["tmax"][0], 30.0)
        self.assertEqual(s["dt_global_c"][0], 0.0)

    def test_contraste_uniforme_robin_coincide_con_formula_cerrada(self):
        grilla, solver = solver_linea([1, 1, 1, 1], robin(8.0, 20.0))
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.full((4, 1, 1), 30.0), solver)
        # Todas las caras tienen la misma T de cara: Rc = 0.25 m, 1/h = 0.125 → T_cara = 30 − 10·(0.25/0.375)
        esperado_cara = temperatura_cara_estimada(robin(8.0, 20.0), 30.0, 2.0, 1.0, 0.0)
        s = extremos.serie("global")
        self.assertAlmostEqual(s["fmin"][0], esperado_cara)
        self.assertAlmostEqual(s["con_signo"][0], 30.0 - esperado_cara)
        self.assertGreater(s["con_signo"][0], 0.0)

    def test_dirichlet_igual_a_la_temperatura_uniforme_da_contraste_cero(self):
        grilla, solver = solver_linea([1, 1, 1, 1], dirichlet(30.0))
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.full((4, 1, 1), 30.0), solver)
        self.assertEqual(extremos.serie("global")["contraste_abs_c"][0], 0.0)


class ExtremosConocidosTest(unittest.TestCase):
    def test_extremos_ubicaciones_y_signo(self):
        grilla, solver = solver_linea([1, 1, 1, 1], dirichlet(10.0))
        extremos = SeguimientoExtremos(grilla)
        campo = np.array([12.0, 40.0, 25.0, -5.0]).reshape(4, 1, 1)
        extremos.registrar(0.0, campo, solver)
        s = extremos.serie("global")
        self.assertEqual((s["tmax"][0], s["tmin"][0]), (40.0, -5.0))
        self.assertEqual(s["dt_global_c"][0], 45.0)
        fila = next(f for f in extremos.resumen() if f["metrica"] == "dT_global_celdas")
        self.assertEqual(fila["punto_a"], (1.5, 0.5, 0.5))      # celda más caliente
        self.assertEqual(fila["punto_b"], (3.5, 0.5, 0.5))      # celda más fría
        # Caras a 10 °C: |40 − 10| = 30 supera a |−5 − 10| = 15; signo T_celda − T_cara > 0.
        self.assertEqual(s["con_signo"][0], 30.0)
        self.assertEqual(s["con_tcel"][0], 40.0)
        self.assertEqual(s["con_tcara"][0], 10.0)
        self.assertEqual(int(s["con_cel"][0]), 1)

    def test_contraste_negativo_cuando_domina_la_celda_fria(self):
        grilla, solver = solver_linea([1, 1, 1, 1], dirichlet(30.0))
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.array([28.0, 31.0, 29.0, -5.0]).reshape(4, 1, 1), solver)
        s = extremos.serie("global")
        self.assertEqual(s["con_signo"][0], -35.0)               # Tcelda_min − Tcara
        self.assertEqual(int(s["con_cel"][0]), 3)

    def test_orientacion_y_posicion_de_la_cara_del_contraste(self):
        grilla, solver = solver_linea([1, 1, 1, 1], neumann(0.0))   # T_cara = T_celda
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.array([12.0, 40.0, 25.0, -5.0]).reshape(4, 1, 1), solver)
        s = extremos.serie("global")
        self.assertEqual(s["con_signo"][0], 45.0)                 # 40 − (−5)
        self.assertEqual(int(s["con_cel"][0]), 1)                 # celda caliente
        self.assertEqual(int(s["con_fcel"][0]), 3)                # cara de la celda fría
        fila = next(f for f in extremos.resumen() if f["metrica"] == "contraste_celda_cara")
        self.assertEqual(fila["cara"], "x+")                      # primera cara expuesta en ese orden
        self.assertEqual(fila["punto_b"], (4.0, 0.5, 0.5))        # centro de la cara, a d/2 de la celda

    def test_contraste_extremos_coincide_con_busqueda_exhaustiva(self):
        rng = np.random.default_rng(3)
        forma = (4, 3, 3)
        ids = np.ones(forma, dtype=np.int32)
        mascara = rng.random(forma) > 0.25
        ids[~mascara] = 0
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (0.5, 0.4, 0.3), ids, {1: 1.0})
        solver = DifusionCartesiana(
            forma=forma, dx_m=(0.5, 0.4, 0.3), conductividad_w_m_k=rng.uniform(1, 3, forma),
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            fronteras={"x-": robin(8.0, 5.0), "x+": dirichlet(15.0), "y-": neumann(30.0),
                       "y+": robin(3.0, 35.0), "z-": neumann(0.0), "z+": robin(12.0, 20.0)},
            mascara_activa=mascara)
        campo = rng.uniform(10, 60, forma)
        campo[~mascara] = 1e6                                     # basura en celdas inactivas
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, campo, solver)
        # Búsqueda exhaustiva sobre TODOS los pares (celda activa, cara expuesta).
        mejor = 0.0
        celdas = list(zip(*np.nonzero(mascara)))
        for indice, cara in itertools.product(celdas, CARAS):
            if not solver.superficie_expuesta[cara][indice]:
                continue
            eje = "xyz".index(cara[0])
            t_cara = temperatura_cara_estimada(solver.fronteras[cara], campo[indice],
                                               solver.k[indice], solver.dx_m[eje], 0.0)
            for otra in celdas:
                mejor = max(mejor, abs(campo[otra] - t_cara))
        s = extremos.serie("global")
        self.assertAlmostEqual(s["contraste_abs_c"][0], mejor, places=9)
        self.assertLess(s["tmax"][0], 1e5)
        self.assertEqual(s["tmin"][0], campo[mascara].min())


class CeldasInactivasTest(unittest.TestCase):
    def test_se_excluyen_celdas_inactivas_de_celdas_y_caras(self):
        grilla, solver = solver_linea([1, 1, 0, 0], dirichlet(10.0))
        campo = np.array([20.0, 30.0, 999.0, -999.0]).reshape(4, 1, 1)
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, campo, solver)
        s = extremos.serie("global")
        self.assertEqual((s["tmin"][0], s["tmax"][0]), (20.0, 30.0))
        self.assertEqual(s["n_celdas"][0], 2)
        # Solo las caras de las 2 celdas activas: x- (1) + y±,z± (2 celdas × 4) + x+ de la celda 1 = 10.
        self.assertEqual(s["n_caras"][0], 10)
        self.assertEqual(s["fmin"][0], 10.0)
        self.assertEqual(s["fmax"][0], 10.0)

    def test_alcance_sin_celdas_activas_es_na(self):
        ids = np.zeros((2, 1, 1), dtype=np.int32)
        ids[0, 0, 0], ids[1, 0, 0] = 1, 2
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 1.0, 2: 1.0})
        solver = DifusionCartesiana(
            forma=(2, 1, 1), dx_m=(1.0,) * 3, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
            calor_especifico_j_kg_k=900.0, fronteras=fronteras(dirichlet(10.0)),
            mascara_activa=np.array([True, False]).reshape(2, 1, 1))
        extremos = SeguimientoExtremos(grilla)
        extremos.registrar(0.0, np.array([20.0, 0.0]).reshape(2, 1, 1), solver)
        s2 = extremos.serie("vertida_2")
        for campo in ("tmin", "tmax", "dt_global_c", "con_signo", "contraste_abs_c"):
            self.assertTrue(np.isnan(s2[campo][0]), campo)       # n/a, no 0
        self.assertFalse(np.isnan(extremos.serie("vertida_1")["tmin"][0]))


class SimulacionTest(unittest.TestCase):
    def setUp(self):
        self.grilla = grilla_dos_capas()
        self.vertidas = [vertida(1, 0.0, 18.0), vertida(2, 6.0, 26.0)]

    def correr(self, **kw):
        extremos = SeguimientoExtremos(self.grilla)
        resultado = simular(self.grilla, self.vertidas, extremos, **kw)
        return extremos, resultado

    def test_vertida_es_na_antes_de_activarse_y_parte_en_su_colocacion(self):
        extremos, _ = self.correr()
        s = extremos.serie("vertida_2")
        horas = s["tiempos_s"] / H
        antes = horas < 6.0
        for campo in ("tmin", "tmax", "dt_global_c", "con_signo", "contraste_abs_c"):
            self.assertTrue(np.all(np.isnan(s[campo][antes])), campo)
        n = int(np.argmax(horas >= 6.0))
        self.assertEqual(horas[n], 6.0)
        self.assertEqual((s["tmin"][n], s["tmax"][n]), (26.0, 26.0))   # recién colocada, uniforme
        self.assertEqual(s["dt_global_c"][n], 0.0)
        self.assertEqual(s["n_celdas"][n], 9)
        # El global ya incluye a la vertida 2 en ese mismo instante (la fila se reemplaza).
        self.assertEqual(extremos.serie("global")["tmax"][n], 26.0)
        self.assertEqual(np.count_nonzero(np.isclose(horas, 6.0)), 1)   # sin instantes duplicados

    def test_cara_superior_de_la_vertida_1_pasa_a_interfaz(self):
        extremos, _ = self.correr()
        horas = extremos.serie("global")["tiempos_s"] / H
        s1, s2, sg = (extremos.serie(a) for a in ("vertida_1", "vertida_2", "global"))
        antes, despues = horas < 6.0, horas >= 6.0
        # Vertida 1: 9 (z-) + 9 (z+) + 12 laterales = 30 caras; con la 2 encima pierde las 9 de arriba.
        self.assertTrue(np.all(s1["n_caras"][antes] == 30))
        self.assertTrue(np.all(s1["n_caras"][despues] == 21))
        self.assertTrue(np.all(s2["n_caras"][despues] == 21))          # pierde su cara inferior
        self.assertTrue(np.all(sg["n_caras"][despues] == 42))
        # Ninguna cara del conjunto de contraste de la vertida 1 es z+ después de activarse la 2.
        despues_idx = np.flatnonzero(despues)
        for n in despues_idx[::40]:
            self.assertNotEqual(int(s1["con_fcara"][n]), CARAS.index("z+"))
            self.assertNotEqual(int(s1["fmin_cara"][n]), CARAS.index("z+"))
            self.assertNotEqual(int(s1["fmax_cara"][n]), CARAS.index("z+"))

    def test_registra_cada_paso_y_cada_activacion(self):
        extremos, _ = self.correr()
        t = np.asarray(extremos.tiempos_s)
        self.assertTrue(np.all(np.diff(t) > 0))
        self.assertLessEqual(np.max(np.diff(t)), 1800.0 + 1e-9)
        self.assertEqual(t[-1], 24 * H)

    def test_no_cambia_la_solucion(self):
        base = simular(self.grilla, self.vertidas)
        extremos = SeguimientoExtremos(self.grilla)
        puntos = SeguimientoPuntos(self.grilla, resolver_puntos(self.grilla, {
            "1": {"interior_m": [1.5, 1.5, 0.5], "voxel_proximo_m": [2.5, 1.5, 0.5],
                  "cara_superficie": "x+"}}))
        con = simular(self.grilla, self.vertidas, [puntos, extremos])      # ambos observadores
        np.testing.assert_array_equal(base.tiempos_s, con.tiempos_s)
        for a, b in zip(base.temperaturas_c, con.temperaturas_c):
            np.testing.assert_array_equal(a, b)
        self.assertEqual(base.energia_final_j, con.energia_final_j)
        self.assertEqual(base.error_balance_j, con.error_balance_j)
        self.assertGreater(len(extremos.tiempos_s), len(base.tiempos_s))   # observa más pasos que los guardados

    def test_maximos_en_todos_los_pasos_superan_o_igualan_a_los_de_fotogramas(self):
        extremos, _ = self.correr()
        marcos = [0.0, 6 * H, 12 * H, 24 * H]
        filas = {(f["alcance"], f["metrica"], f["base"]): f for f in extremos.resumen(marcos)}
        for alcance in extremos.alcances:
            for metrica in ("dT_global_celdas", "contraste_celda_cara"):
                completo = filas[(alcance, metrica, "todos_los_pasos")]
                marco = filas[(alcance, metrica, "fotogramas")]
                if marco["valido"]:
                    self.assertLessEqual(marco["abs_c"], completo["abs_c"] + 1e-12)
                    self.assertTrue(any(math.isclose(marco["tiempo_h"], m / H, abs_tol=1e-6) for m in marcos))
        # En esta corrida el máximo de contraste global ocurre entre fotogramas.
        glo = filas[("global", "contraste_celda_cara", "todos_los_pasos")]
        self.assertFalse(any(math.isclose(glo["tiempo_h"], m / H, abs_tol=1e-6) for m in marcos))
        self.assertGreater(glo["abs_c"], filas[("global", "contraste_celda_cara", "fotogramas")]["abs_c"])

    def test_csv_usa_na_antes_de_la_activacion(self):
        extremos, _ = self.correr()
        with tempfile.TemporaryDirectory() as carpeta:
            ruta = Path(carpeta) / "evolucion.csv"
            extremos.guardar_evolucion_csv(ruta)
            extremos.guardar_maximos_csv(Path(carpeta) / "maximos.csv", [0.0, 24 * H])
            with ruta.open(encoding="utf-8") as archivo:
                filas = list(csv.DictReader(archivo))
            with (Path(carpeta) / "maximos.csv").open(encoding="utf-8") as archivo:
                maximos = list(csv.DictReader(archivo))
        previas = [f for f in filas if f["alcance"] == "vertida_2" and float(f["tiempo_h"]) < 6.0]
        self.assertTrue(previas)
        for fila in previas:
            self.assertEqual(fila["Tmin_celdas_C"], "n/a")
            self.assertEqual(fila["dT_global_celdas_C"], "n/a")
            self.assertEqual(fila["contraste_abs_celda_cara_C"], "n/a")
        self.assertEqual({m["base_de_busqueda"] for m in maximos}, {"todos_los_pasos", "fotogramas"})
        self.assertEqual(len(maximos), 3 * 2 * 2)             # 3 alcances × 2 métricas × 2 bases


if __name__ == "__main__":
    unittest.main()
