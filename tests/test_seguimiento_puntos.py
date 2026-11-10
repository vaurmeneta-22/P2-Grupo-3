"""Pruebas del seguimiento de temperaturas y diferencias por vertida."""
import math
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.diffusion_nd import dirichlet, neumann, robin
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.seguimiento import (
    SeguimientoPuntos,
    proponer_puntos,
    resolver_puntos,
    resumen_diferencias,
    temperatura_cara_estimada,
)
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
    ids[:, :, 0] = 1
    ids[:, :, 1] = 2
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 9.0, 2: 9.0})


def simular(grilla, vertidas, seguimiento=None, horas=24.0):
    return simular_vertidas(
        grilla=grilla, vertidas=vertidas, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
        calor_especifico_j_kg_k=900.0, hidratacion=ModeloPulsoCosenoidal(500.0, 12 * H),
        tiempo_final_s=horas * H, paso_tiempo_max_s=1800.0, tiempos_muestreo_s=[8 * H],
        fronteras={cara: robin(8.0, 20.0) for cara in ("x-", "x+", "y-", "y+", "z-", "z+")},
        seguimiento=seguimiento)


PUNTOS = {
    "1": {"interior_m": [1.5, 1.5, 0.5], "voxel_proximo_m": [2.5, 1.5, 0.5], "cara_superficie": "x+"},
    "2": {"interior_m": [1.5, 1.5, 1.5], "voxel_proximo_m": [2.5, 1.5, 1.5], "cara_superficie": "x+"},
}
PUNTOS_TOPE = {  # cara superior de la vertida 1: pasa a interfaz cuando se activa la 2
    "1": {"interior_m": [1.5, 1.5, 0.5], "voxel_proximo_m": [1.5, 1.5, 0.5], "cara_superficie": "z+"},
}


class TemperaturaCaraTest(unittest.TestCase):
    def test_dirichlet_es_la_prescrita(self):
        self.assertEqual(temperatura_cara_estimada(dirichlet(30.0), 55.0, 2.0, 0.5, 0.0), 30.0)

    def test_neumann_usa_flujo_entrante_y_media_celda(self):
        # q entrante calienta la cara respecto del centro: T_cara = T + q d/(2k)
        valor = temperatura_cara_estimada(neumann(100.0), 20.0, 2.0, 0.5, 0.0)
        self.assertAlmostEqual(valor, 20.0 + 100.0 * 0.5 / 4.0)
        self.assertEqual(temperatura_cara_estimada(neumann(0.0), 20.0, 2.0, 0.5, 0.0), 20.0)

    def test_robin_perfil_lineal_estacionario(self):
        """Barra 1D en régimen: la cara reconstruida coincide con la extrapolación lineal."""
        k, h, t_amb, t_izq, n, d = 2.0, 5.0, 10.0, 50.0, 10, 0.3
        largo = n * d
        q = (t_izq - t_amb) / (largo / k + 1.0 / h)            # flujo hacia +x [W/m²]
        x_centros = (np.arange(n) + 0.5) * d
        perfil = t_izq - q * x_centros / k                     # Dirichlet 50 °C en x=0
        esperado = t_izq - q * largo / k                       # T(L)
        obtenido = temperatura_cara_estimada(robin(h, t_amb), perfil[-1], k, d, 0.0)
        self.assertAlmostEqual(obtenido, esperado, places=10)
        # La cara debe quedar entre la celda y el ambiente.
        self.assertTrue(t_amb < obtenido < perfil[-1])

    def test_robin_coincide_con_el_flujo_del_solver(self):
        """T_cara - T_celda = q_in Rc, con q_in = G_superficie (T_amb - T)/A del solver."""
        from mcplan.thermal.diffusion_nd import DifusionCartesiana
        solver = DifusionCartesiana(
            forma=(2, 2, 2), dx_m=(0.5, 0.5, 0.5), conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            fronteras={"x+": robin(8.0, 20.0)}, mascara_activa=np.ones((2, 2, 2), dtype=bool))
        t_celda = 33.0
        q_in = solver.G_superficie["x+"][1, 0, 0] / solver.area_cara_m2[0] * (20.0 - t_celda)
        esperado = t_celda + q_in * 0.5 / (2.0 * 2.0)
        self.assertAlmostEqual(
            temperatura_cara_estimada(robin(8.0, 20.0), t_celda, 2.0, 0.5, 0.0), esperado)


class ResumenDiferenciasTest(unittest.TestCase):
    def test_maximo_absoluto_conserva_signo_e_instante(self):
        tiempos = np.array([0.0, 1.0, 2.0, 3.0, 4.0]) * H
        diferencia = np.array([math.nan, 2.0, -5.0, 4.0, math.nan])
        resumen = resumen_diferencias(tiempos, diferencia)
        self.assertEqual(resumen["max_abs_c"], 5.0)
        self.assertEqual(resumen["diferencia_con_signo_c"], -5.0)
        self.assertEqual(resumen["tiempo_h"], 2.0)

    def test_sin_muestras_validas_es_none(self):
        self.assertIsNone(resumen_diferencias([0.0, 1.0], [math.nan, math.nan]))


class ResolverPuntosTest(unittest.TestCase):
    def test_punto_en_otra_vertida_falla(self):
        malos = {"1": dict(PUNTOS["1"], interior_m=[1.5, 1.5, 1.5])}   # cae en la vertida 2
        with self.assertRaisesRegex(ValueError, "pertenece a VertidaID=2"):
            resolver_puntos(grilla_dos_capas(), malos)

    def test_punto_fuera_de_la_grilla_falla(self):
        malos = {"1": dict(PUNTOS["1"], voxel_proximo_m=[9.5, 1.5, 0.5])}
        with self.assertRaisesRegex(ValueError, "fuera de la grilla"):
            resolver_puntos(grilla_dos_capas(), malos)

    def test_cara_interna_a_la_vertida_falla(self):
        malos = {"1": dict(PUNTOS["1"], voxel_proximo_m=[1.5, 1.5, 0.5], cara_superficie="x+")}
        with self.assertRaisesRegex(ValueError, "interna"):
            resolver_puntos(grilla_dos_capas(), malos)

    def test_informa_vecino_y_geometria(self):
        punto = resolver_puntos(grilla_dos_capas(), PUNTOS_TOPE)[0]
        self.assertEqual(punto.vecino_id, 2)
        self.assertEqual(punto.centro_cara_m, (1.5, 1.5, 1.0))
        self.assertAlmostEqual(punto.distancia_interior_a_cara_m, 0.5)

    def test_propuesta_es_valida_y_pertenece_a_cada_vertida(self):
        grilla = grilla_dos_capas()
        propuesta = proponer_puntos(grilla)
        for punto in resolver_puntos(grilla, propuesta):
            self.assertEqual(grilla.vertida_id[punto.indice_interior], punto.vertida_id)
            self.assertEqual(grilla.vertida_id[punto.indice_voxel_proximo], punto.vertida_id)
            self.assertEqual(punto.vecino_id, 0)


def grilla_bloques(refinamiento):
    """Dos vertidas apiladas de 8x8 m en planta (z 0-2 y 2-4 m), con dx = 1/refinamiento."""
    n = refinamiento
    ids = np.zeros((8 * n, 8 * n, 4 * n), dtype=np.int32)
    ids[:, :, :2 * n] = 1
    ids[:, :, 2 * n:] = 2
    paso = 1.0 / n
    return GrillaVoxel((0.0, 0.0, 0.0), (paso,) * 3, ids, {1: 128.0, 2: 128.0})


class ModoAutomaticoTest(unittest.TestCase):
    def test_propuesta_valida_para_distintos_pasos_de_voxel(self):
        for refinamiento in (1, 2, 4):
            grilla = grilla_bloques(refinamiento)
            puntos = resolver_puntos(grilla, proponer_puntos(grilla))   # misma validación estricta
            self.assertEqual([p.vertida_id for p in puntos], [1, 2])
            for p in puntos:
                self.assertEqual(p.cara, "x+")
                self.assertEqual(p.vecino_id, 0)                        # cara libre
                self.assertEqual(grilla.vertida_id[p.indice_interior], p.vertida_id)
                self.assertEqual(grilla.vertida_id[p.indice_voxel_proximo], p.vertida_id)
                self.assertAlmostEqual(p.centro_cara_m[0], 8.0)         # cara exterior del bloque

    def test_las_posiciones_cambian_con_la_resolucion(self):
        centros = [resolver_puntos(grilla_bloques(r), proponer_puntos(grilla_bloques(r)))[0]
                   .centro_interior_m for r in (1, 2)]
        self.assertNotEqual(centros[0], centros[1])      # no es el mismo punto físico

    def test_cara_none_elige_una_lateral_libre(self):
        grilla = grilla_dos_capas()
        for punto in resolver_puntos(grilla, proponer_puntos(grilla, None)).__iter__():
            self.assertIn(punto.cara, ("x+", "x-", "y+", "y-"))
            self.assertEqual(punto.vecino_id, 0)

    def test_cara_forzada_ocupada_falla(self):
        ids = np.zeros((4, 1, 1), dtype=np.int32)
        ids[:2, 0, 0], ids[2:, 0, 0] = 1, 2           # la vertida 1 tiene a la 2 contra x+
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 2.0, 2: 2.0})
        with self.assertRaisesRegex(ValueError, "no se encontró una cara lateral libre"):
            proponer_puntos(grilla, "x+")
        self.assertEqual(proponer_puntos(grilla, None)["1"]["cara_superficie"], "x-")

    def test_csv_separa_centro_del_voxel_y_centro_de_la_cara(self):
        import csv
        import tempfile
        from pathlib import Path
        grilla = grilla_dos_capas()
        seguimiento = SeguimientoPuntos(grilla, resolver_puntos(grilla, PUNTOS))
        simular(grilla, [vertida(1, 0.0, 18.0), vertida(2, 6.0, 26.0)], seguimiento)
        with tempfile.TemporaryDirectory() as carpeta:
            ruta = Path(carpeta) / "s.csv"
            seguimiento.guardar_csv(ruta)
            with ruta.open(encoding="utf-8") as archivo:
                fila = next(csv.DictReader(archivo))
        self.assertEqual(float(fila["x_centro_voxel_proximo_m"]), 2.5)
        self.assertEqual(float(fila["x_centro_cara_T_superficie_m"]), 3.0)   # a d/2 = 0.5 m
        self.assertEqual(float(fila["x_centro_voxel_interior_m"]), 1.5)


class SeguimientoSimulacionTest(unittest.TestCase):
    def setUp(self):
        self.grilla = grilla_dos_capas()
        self.vertidas = [vertida(1, 0.0, 18.0), vertida(2, 6.0, 26.0)]

    def seguir(self, puntos):
        seguimiento = SeguimientoPuntos(self.grilla, resolver_puntos(self.grilla, puntos))
        resultado = simular(self.grilla, self.vertidas, seguimiento)
        return seguimiento, resultado

    def test_antes_de_activarse_es_na_y_al_activarse_parte_en_colocacion(self):
        seguimiento, _ = self.seguir(PUNTOS)
        serie = seguimiento.serie(2)
        horas = serie["tiempos_s"] / H
        antes = horas < 6.0
        self.assertTrue(np.all(np.isnan(serie["interior_c"][antes])))
        self.assertTrue(np.all(np.isnan(serie["dt_cara_c"][antes])))
        self.assertTrue(all(e.startswith("n/a") for e in np.asarray(serie["estado"])[antes]))
        primero = int(np.argmax(horas >= 6.0))
        self.assertEqual(horas[primero], 6.0)
        self.assertEqual(serie["interior_c"][primero], 26.0)       # no 0 °C ni otra temperatura
        self.assertEqual(serie["celda_proxima_c"][primero], 26.0)
        self.assertEqual(serie["estado"][primero], "expuesta")
        # La vertida 1 ya estaba activa en t=0 con su temperatura de colocación.
        self.assertEqual(seguimiento.serie(1)["interior_c"][0], 18.0)

    def test_registra_cada_paso_de_integracion(self):
        seguimiento, _ = self.seguir(PUNTOS)
        tiempos = np.asarray(seguimiento.tiempos_s)
        self.assertTrue(np.all(np.diff(tiempos) > 0))              # sin instantes duplicados
        self.assertLessEqual(np.max(np.diff(tiempos)), 1800.0 + 1e-9)
        self.assertEqual(tiempos[-1], 24 * H)

    def test_diferencia_es_interior_menos_superficie(self):
        seguimiento, _ = self.seguir(PUNTOS)
        serie = seguimiento.serie(1)
        validos = np.isfinite(serie["dt_cara_c"])
        np.testing.assert_allclose(serie["dt_cara_c"][validos],
                                   (serie["interior_c"] - serie["cara_estimada_c"])[validos])
        np.testing.assert_allclose(serie["dt_celda_c"][validos],
                                   (serie["interior_c"] - serie["celda_proxima_c"])[validos])
        fila = seguimiento.resumen()[0]["dt_cara"]
        self.assertAlmostEqual(fila["max_abs_c"], float(np.nanmax(np.abs(serie["dt_cara_c"]))))

    def test_cara_pasa_a_interfaz_al_activarse_la_vecina(self):
        seguimiento, _ = self.seguir(PUNTOS_TOPE)
        serie = seguimiento.serie(1)
        horas = serie["tiempos_s"] / H
        estados = np.asarray(serie["estado"])
        self.assertTrue(np.all(estados[horas < 6.0] == "expuesta (vecina futura: vertida 2)"))
        despues = horas >= 6.0
        self.assertTrue(np.all(estados[despues] == "interfaz con vertida 2 (activa)"))
        self.assertTrue(np.all(np.isnan(serie["cara_estimada_c"][despues])))   # ya no hay superficie
        self.assertTrue(np.all(np.isnan(serie["dt_cara_c"][despues])))
        self.assertTrue(np.all(np.isfinite(serie["celda_proxima_c"])))          # la celda sigue activa
        cambios = seguimiento.cambios_de_estado(1)
        self.assertEqual([round(h, 6) for h, _ in cambios], [0.0, 6.0])

    def test_el_seguimiento_no_cambia_la_solucion(self):
        base = simular(self.grilla, self.vertidas)
        con_seguimiento, _ = self.seguir(PUNTOS)
        _, resultado = self.seguir(PUNTOS)
        self.assertIsNotNone(con_seguimiento)
        np.testing.assert_array_equal(base.tiempos_s, resultado.tiempos_s)
        for campo_base, campo in zip(base.temperaturas_c, resultado.temperaturas_c):
            np.testing.assert_array_equal(campo_base, campo)
        np.testing.assert_array_equal(base.temperaturas_maximas_c, resultado.temperaturas_maximas_c)
        self.assertEqual(base.energia_final_j, resultado.energia_final_j)
        self.assertEqual(base.error_balance_j, resultado.error_balance_j)

    def test_csv_usa_na_y_no_cero(self):
        import csv
        import tempfile
        from pathlib import Path
        seguimiento, _ = self.seguir(PUNTOS)
        with tempfile.TemporaryDirectory() as carpeta:
            ruta = Path(carpeta) / "seguimiento.csv"
            seguimiento.guardar_csv(ruta)
            with ruta.open(encoding="utf-8") as archivo:
                filas = list(csv.DictReader(archivo))
            seguimiento.guardar_resumen_csv(Path(carpeta) / "resumen.csv")
        previas = [f for f in filas if f["VertidaID"] == "2" and float(f["tiempo_h"]) < 6.0]
        self.assertTrue(previas)
        for fila in previas:
            self.assertEqual(fila["T_interior_C"], "n/a")
            self.assertEqual(fila["dT_cara_C"], "n/a")


if __name__ == "__main__":
    unittest.main()
