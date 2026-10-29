"""Pruebas locales del ejercicio 3D de fuente variable (solo usan el solver_ejercicio.py de esta carpeta).

Ejecutar:  python -m unittest discover -s ejercicios_locales/fuente_variable_3d -v
Usan mallas gruesas para ser rapidas. No forman parte de la bateria de la entrega (tests/).
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import csv
import math
import os
import sys
import unittest

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np

CARPETA = Path(__file__).resolve().parent


def _cargar(nombre, ruta):
    spec = spec_from_file_location(nombre, ruta)
    modulo = module_from_spec(spec)
    sys.modules[nombre] = modulo
    spec.loader.exec_module(modulo)
    return modulo


ej = _cargar("ejercicio_fuente_variable_3d_mod", CARPETA / "ejercicio_fuente_variable_3d.py")
solver = ej.solver
GRUESA = (12, 12, 3)


class TestCoherenciaConElEjercicio2D(unittest.TestCase):
    def test_propiedades_igual_al_2d(self):
        ruta = CARPETA.parent / "fuente_variable_2d" / "ejercicio_fuente_variable_2d.py"
        if not ruta.exists():
            self.skipTest("no esta el ejercicio 2D")
        d2 = _cargar("ejercicio_fuente_variable_2d_ref3", ruta)
        for n in ("CONDUCTIVIDAD_W_M_K", "DENSIDAD_KG_M3", "CALOR_ESPECIFICO_J_KG_K", "TEMPERATURA_INICIAL_C", "DURACION_FUENTE_S"):
            self.assertEqual(getattr(ej, n), getattr(d2, n), n)
        for t in (0.0, 1234.5, 86400.0, 172800.0, 200000.0):
            self.assertEqual(ej.fuente_q(t, 55.0), d2.fuente_q(t, 55.0))

    def test_qmax_coincide_con_la_corrida_dirichlet_2d_si_existe(self):
        ruta = CARPETA.parent / "fuente_variable_2d" / "resultados" / "centro_cada_paso.csv"
        if not ruta.exists():
            self.skipTest("no hay resultados del ejercicio 2D (se regeneran ejecutandolo)")
        with open(ruta, encoding="utf-8") as f:
            filas = list(csv.DictReader(f))
        est = []
        for fila in filas:
            t, q = float(fila["tiempo_s"]), float(fila["Q_inicio_de_paso_W_m3"])
            g = 0.5 * (1.0 - math.cos(2.0 * math.pi * t / ej.DURACION_FUENTE_S))
            if 0.3 < g <= 1.0:
                est.append(q / g)
        self.assertAlmostEqual(float(np.mean(est)), ej.QMAX_W_M3, delta=1e-9)


class TestConfiguracion(unittest.TestCase):
    def test_h_igual_k_beta_y_nombres(self):
        self.assertEqual(ej.h_desde_beta(1.0), ej.CONDUCTIVIDAD_W_M_K)
        self.assertEqual(ej.h_desde_beta(0.1, 5.0), 0.5)
        self.assertEqual(ej.nombre_beta(0.1), "robin_beta_0p1")
        esc = {e["slug"]: e for e in ej.definir_escenarios()}
        self.assertEqual(set(esc), {"dirichlet_20C", "dirichlet_30C", "neumann_cero", "robin_beta_0p1", "robin_beta_1", "robin_beta_10"})
        c = esc["robin_beta_10"]["fronteras"]["x-"]
        self.assertEqual((c.tipo, c.coeficiente_h_w_m2_k, c.temperatura(0.0)), ("robin", 20.0, 20.0))
        self.assertEqual(len(esc["dirichlet_30C"]["fronteras"]), 6)

    def test_dt_comun_es_estable_en_todos_los_escenarios_con_aristas_y_esquinas(self):
        escenarios = ej.definir_escenarios()
        dt, pasos, limites = ej.paso_comun(escenarios)
        self.assertLessEqual(dt, ej.FRACCION_PASO_ESTABLE * min(limites.values()) * (1 + 1e-9))
        # Dirichlet en 6 caras: la esquina (3 caras Dirichlet + 3 internas) manda: rho cp dx² / (9 k)
        dx = ej.LADO_X_M / ej.CELDAS_X
        self.assertAlmostEqual(limites["dirichlet_20C"], ej.RHOCP * dx**2 / (9 * ej.CONDUCTIVIDAD_W_M_K), places=6)
        self.assertGreater(limites["neumann_cero"], limites["dirichlet_20C"])
        self.assertAlmostEqual(ej.TIEMPO_FINAL_S / (ej.NUM_FOTOGRAMAS - 1) / pasos, dt)


class TestCortes(unittest.TestCase):
    def test_corte_interpolado_exacto_en_campo_lineal(self):
        celdas, lados = (6, 6, 4), (3.0, 3.0, 0.6)
        dx = tuple(L / n for L, n in zip(lados, celdas))
        ejes = [(np.arange(n) + 0.5) * d for n, d in zip(celdas, dx)]
        X, Y, Z = np.meshgrid(*ejes, indexing="ij")
        T = 1.0 + 2.0 * X + 3.0 * Y + 5.0 * Z
        coef = (2.0, 3.0, 5.0)
        for eje, coord in ((2, 0.30), (1, 1.5), (0, 1.5), (2, 0.2)):
            plano, desc = ej.corte_interpolado(T, eje, coord, dx)
            resto = [i for i in range(3) if i != eje]
            A, B = np.meshgrid(ejes[resto[0]], ejes[resto[1]], indexing="ij")
            esperado = 1.0 + coef[resto[0]] * A + coef[resto[1]] * B + coef[eje] * coord
            np.testing.assert_allclose(plano, esperado, rtol=1e-12)
            self.assertIn(f"{'xyz'[eje]} = {coord:.3f} m", desc)
        self.assertIn("interpolado", ej.corte_interpolado(T, 2, 0.30, dx)[1])      # plano 0.30 m cae entre capas 0.225 y 0.375


class TestVerificaciones(unittest.TestCase):
    def test_verificaciones_rapidas_del_script(self):
        for nombre, ok, detalle in ej.verificaciones_rapidas():
            self.assertTrue(ok, f"{nombre}: {detalle}")

    def test_neumann_cero_uniforme_71_30_y_balance(self):
        esc = {e["slug"]: e for e in ej.definir_escenarios()}["neumann_cero"]
        dt = 0.8 * ej.construir_sistema(esc["fronteras"], GRUESA).paso_maximo_estable_s()
        s, r, _ = ej.simular(esc, dt, GRUESA)
        self.assertLess(float(r.temperatura_c.max() - r.temperatura_c.min()), 1e-9)
        self.assertAlmostEqual(float(r.temperatura_c.mean()), 71.30, delta=0.01)
        self.assertEqual(r.energia_frontera_total_j, 0.0)
        self.assertAlmostEqual(r.energia_final_j - r.energia_inicial_j, r.energia_fuente_j, delta=1e-9 * r.energia_fuente_j)
        err = float(np.max(np.abs(r.historial_puntos_c[:, 0] - ej.temperatura_adiabatica_exacta(r.historial_tiempos_s))))
        self.assertLessEqual(err, ej.tolerancia_integracion_c(r.paso_tiempo_s))

    def test_simetria_metricas_y_balance_robin(self):
        esc = {e["slug"]: e for e in ej.definir_escenarios()}["robin_beta_1"]
        dt = 0.8 * ej.construir_sistema(esc["fronteras"], GRUESA).paso_maximo_estable_s()
        s, r, _ = ej.simular(esc, dt, GRUESA)
        self.assertLess(ej.simetria_max_dif(r.temperatura_c), 1e-9)
        self.assertEqual(len(r.historial_max_c), r.numero_pasos + 1)
        self.assertTrue(np.all(r.historial_max_c >= r.historial_puntos_c[:, 0] - 1e-9))     # el maximo global acota al central
        self.assertTrue(np.all(r.historial_min_c <= r.historial_media_c + 1e-12))
        self.assertLess(abs(r.error_balance_j), 1e-9 * abs(r.energia_fuente_j))
        self.assertLess(r.energia_frontera_total_j, r.energia_fuente_j)

    def test_dirichlet_30_calienta_y_dirichlet_20_enfria_respecto_del_aislado(self):
        esc = {e["slug"]: e for e in ej.definir_escenarios()}
        dt = 0.8 * min(ej.construir_sistema(e["fronteras"], GRUESA).paso_maximo_estable_s() for e in esc.values())
        tc = {k: ej.simular(esc[k], dt, GRUESA)[1].historial_puntos_c[:, 0].max() for k in ("dirichlet_20C", "dirichlet_30C", "neumann_cero")}
        self.assertLess(tc["dirichlet_20C"], tc["dirichlet_30C"])
        self.assertLess(tc["dirichlet_30C"], tc["neumann_cero"])


if __name__ == "__main__":
    unittest.main()
