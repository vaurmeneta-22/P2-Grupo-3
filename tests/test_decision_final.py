"""Pruebas de curvas adiabáticas, S25, planes y agregación probabilística."""
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.decision import (
    PlanConstruccion, ResultadoEscenario, SeguimientoMetricasDecision,
    agregar_probabilidades,
)
from mcplan.thermal.diffusion_nd import DifusionCartesiana, dirichlet, neumann, robin
from mcplan.thermal.hidratacion import ModeloCurvaAdabatica, ModeloPulsoCosenoidal
from mcplan.thermal.vertidas import simular_vertidas


def tetra(identificador=1):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    triangulos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", "gid", 0, 18,
                      "M1", vertices, triangulos, {})


class CurvaAdabaticaTest(unittest.TestCase):
    def test_derivada_por_tramos_y_cola_cero(self):
        curva = ModeloCurvaAdabatica((0.0, 3600.0, 7200.0), (0.0, 2.0, 2.0))
        rho, cp = 2400.0, 900.0
        q = curva.generacion_w_m3(np.array([-1.0, 0.0, 3599.0, 3600.0, 7200.0]), rho, cp)
        self.assertAlmostEqual(q[0], 0.0)
        self.assertAlmostEqual(q[1], rho * cp * 2.0 / 3600.0)
        self.assertAlmostEqual(q[2], rho * cp * 2.0 / 3600.0)
        self.assertEqual(q[3], 0.0)
        self.assertEqual(q[4], 0.0)

    def test_rechaza_curva_no_monotona_o_tiempo_repetido(self):
        with self.assertRaises(ValueError):
            ModeloCurvaAdabatica((0.0, 1.0), (2.0, 1.0))
        with self.assertRaises(ValueError):
            ModeloCurvaAdabatica((0.0, 0.0), (0.0, 1.0))

    def test_ajuste_isotonico_corrige_ruido_sin_modificar_la_curva_cruda(self):
        ajustada = ModeloCurvaAdabatica._ajuste_isotonico((0.0, 2.0, 1.0, 3.0))
        self.assertEqual(ajustada, (0.0, 1.5, 1.5, 3.0))

    def test_curva_adiabatica_funciona_directamente_como_fuente_del_solver(self):
        forma = (4, 1, 1)
        ids = np.ones(forma, dtype=np.int32)
        grilla = GrillaVoxel((0, 0, 0), (0.1, 0.1, 0.1), ids, {1: 0.004})
        resultado = simular_vertidas(
            grilla=grilla, vertidas=[tetra()], conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            hidratacion=ModeloCurvaAdabatica((0.0, 3600.0, 7200.0), (0.0, 1.0, 1.0)),
            tiempo_final_s=3600.0, paso_tiempo_max_s=600.0,
            fronteras={cara: neumann(0.0) for cara in ("x-", "x+", "y-", "y+", "z-", "z+")})
        self.assertAlmostEqual(float(resultado.temperaturas_c[-1].mean()), 19.0, places=6)
        self.assertAlmostEqual(resultado.error_balance_j, 0.0, places=6)


class S25Test(unittest.TestCase):
    def solver(self, dx=0.1):
        forma = (4, 1, 1)
        ids = np.ones(forma, dtype=np.int32)
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (dx, dx, dx), ids, {1: 4 * dx**3})
        solver = DifusionCartesiana(
            forma=forma, dx_m=(dx, dx, dx), conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            fronteras={cara: dirichlet(0.0) for cara in ("x-", "x+", "y-", "y+", "z-", "z+")},
            mascara_activa=np.ones(forma, dtype=bool))
        return grilla, solver

    def test_interpolacion_25mm_y_extremos_globales(self):
        grilla, solver = self.solver()
        observador = SeguimientoMetricasDecision(grilla, 0.025)
        campo = np.array([20.0, 30.0, 40.0, 50.0]).reshape(4, 1, 1)
        observador.registrar(0.0, campo, solver)
        r = observador.resumen()
        # La cara x- vale 0 °C y el centro del voxel está a 50 mm: a 25 mm, T=10 °C.
        self.assertEqual(r["M_T_c"], 50.0)
        self.assertAlmostEqual(r["M_delta_c"], 40.0)
        self.assertEqual(r["t_M_delta_h"], 0.0)
        self.assertAlmostEqual(r["xyz_S25_min_m"][0], 0.025)

    def test_actualiza_superficies_s25_al_activar_otra_vertida(self):
        dx = 0.1
        forma = (2, 1, 1)
        ids = np.array([1, 2], dtype=np.int32).reshape(forma)
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (dx, dx, dx), ids,
                             {1: dx**3, 2: dx**3})
        solver = DifusionCartesiana(
            forma=forma, dx_m=(dx, dx, dx), conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            fronteras={cara: dirichlet(0.0) for cara in
                       ("x-", "x+", "y-", "y+", "z-", "z+")},
            mascara_activa=np.array([True, False]).reshape(forma))
        observador = SeguimientoMetricasDecision(grilla, 0.025)
        campo = np.array([30.0, 20.0]).reshape(forma)

        observador.registrar(0.0, campo, solver)
        solver.actualizar_mascara_activa(np.ones(forma, dtype=bool))
        observador.registrar(3600.0, campo, solver)

        resultado = observador.resumen()
        # La cara x+ de V1 pasa a ser interfaz; la nueva cara x+ de V2
        # queda expuesta y es la S25 más fría (T=10 °C).
        self.assertAlmostEqual(resultado["M_delta_c"], 20.0)
        self.assertAlmostEqual(resultado["t_M_delta_h"], 1.0)
        self.assertAlmostEqual(resultado["xyz_S25_min_m"][0], 0.175)

    def test_actualiza_superficies_s25_al_instalar_o_retirar_proteccion(self):
        dx = 0.1
        forma = (2, 1, 1)
        ids = np.ones(forma, dtype=np.int32)
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (dx, dx, dx), ids, {1: 2 * dx**3})
        caras = ("x-", "x+", "y-", "y+", "z-", "z+")
        proteccion = {cara: np.full(forma, 0.1) for cara in caras}
        solver = DifusionCartesiana(
            forma=forma, dx_m=(dx, dx, dx), conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            fronteras={cara: robin(10.0, 0.0) for cara in caras},
            mascara_activa=np.ones(forma, dtype=bool),
            resistencias_borde=proteccion)
        observador = SeguimientoMetricasDecision(grilla, 0.025)
        campo = np.full(forma, 20.0)

        observador.registrar(0.0, campo, solver)
        solver.actualizar_bordes(resistencias_borde={})
        observador.registrar(3600.0, campo, solver)

        resultado = observador.resumen()
        # Al retirar la protección, la cara reconstruida se enfría y S25
        # debe recalcularse con la nueva resistencia (R_extra = 0).
        self.assertAlmostEqual(resultado["M_delta_c"], 2.0)
        self.assertAlmostEqual(resultado["t_M_delta_h"], 1.0)

    def test_rechaza_malla_que_no_resuelve_25mm(self):
        grilla, solver = self.solver(dx=0.04)
        with self.assertRaisesRegex(ValueError, "no resuelve S25"):
            SeguimientoMetricasDecision(grilla, 0.025).registrar(0.0, np.ones(grilla.forma), solver)

    def test_observador_mide_todos_los_pasos_aceptados(self):
        grilla, _ = self.solver(dx=0.1)
        observer = SeguimientoMetricasDecision(grilla, 0.025)
        simular_vertidas(
            grilla=grilla, vertidas=[tetra()], conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            hidratacion=ModeloPulsoCosenoidal(0.0, 3600.0), tiempo_final_s=30.0,
            paso_tiempo_max_s=10.0, fronteras={c: dirichlet(18.0) for c in
                ("x-", "x+", "y-", "y+", "z-", "z+")},
            temperatura_inicial_c=30.0, seguimiento=observer)
        self.assertGreater(observer.resumen()["n_muestras"], 2)


class PlanYProbabilidadesTest(unittest.TestCase):
    def test_plan_aplica_horas_y_temperaturas_sin_mutar_ifc(self):
        original = [tetra(1), tetra(2)]
        plan = PlanConstruccion("propuesto", {1: 2.0, 2: 8.0}, {1: 3.0, 2: 4.0},
                                {2: 22.0}, frentes_simultaneos_max=1)
        aplicado = plan.aplicar(original)
        self.assertEqual([v.inicio_h for v in aplicado], [2.0, 8.0])
        self.assertEqual(aplicado[1].temperatura_colocacion_c, 22.0)
        self.assertEqual(original[1].inicio_h, 0.0)
        self.assertEqual(plan.duracion_obra_h(t0_h=1.0), 11.0)

    def test_union_no_cuenta_dos_veces_escenario_que_supera_ambos_limites(self):
        filas = [
            ResultadoEscenario("p", "a", 0.04, 71.0, 21.0),
            ResultadoEscenario("p", "b", 0.96, 60.0, 10.0),
        ]
        r = agregar_probabilidades(filas)[0]
        self.assertAlmostEqual(r["P_T"], 0.04)
        self.assertAlmostEqual(r["P_delta"], 0.04)
        self.assertAlmostEqual(r["P_f"], 0.04)
        self.assertTrue(r["cumple_P_f"])

    def test_exige_pesos_unitarios_y_escenarios_unicos(self):
        with self.assertRaisesRegex(ValueError, "deben sumar 1"):
            agregar_probabilidades([ResultadoEscenario("p", "a", 0.5, 10.0, 5.0)])
        with self.assertRaisesRegex(ValueError, "duplicados"):
            agregar_probabilidades([
                ResultadoEscenario("p", "a", 0.5, 10.0, 5.0),
                ResultadoEscenario("p", "a", 0.5, 10.0, 5.0),
            ])


if __name__ == "__main__":
    unittest.main()
