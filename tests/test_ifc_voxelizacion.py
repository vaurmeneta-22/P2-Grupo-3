"""Pruebas iniciales de extracción IFC y voxelización 3D."""
from pathlib import Path
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC, extraer_vertidas
from mcplan.geometry.voxelizacion import (
    area_interfaces_entre_vertidas,
    caras_superficie_exterior,
    voxelizar_vertidas,
)


RAIZ = Path(__file__).resolve().parents[1]
IFC_BASE = RAIZ / "models" / "ifc" / "Fundacion_1_con_propiedades.ifc"


def vertida_cubo(lado_m=1.0):
    """Malla cerrada triangulada de un cubo [0,lado]^3 para prueba unitaria."""
    v = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float) * lado_m
    f = np.array([
        [0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
        [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
        [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7],
    ], dtype=np.int64)
    return VertidaIFC(1, "Vertida_01", "cube-guid", 0.0, 18.0, "M1", v, f,
                      {"VertidaID": 1, "Inicio_h": 0.0,
                       "TemperaturaColocacion_C": 18.0, "TipoHormigon": "M1"})


class ExtraccionIfcTests(unittest.TestCase):
    def test_extrae_las_cinco_vertidas_y_sus_datos(self):
        vertidas = extraer_vertidas(IFC_BASE)
        self.assertEqual([v.vertida_id for v in vertidas], [1, 2, 3, 4, 5])
        self.assertEqual([v.inicio_h for v in vertidas], [0, 6, 12, 18, 24])
        self.assertTrue(all(v.temperatura_colocacion_c == 18.0 for v in vertidas))
        np.testing.assert_allclose(vertidas[-1].dimensiones_m, [7, 7, 1], atol=1e-5)
        self.assertAlmostEqual(max(v.maximo_m[2] for v in vertidas), 4.3, places=5)

    def test_rechaza_un_ifc_inexistente(self):
        with self.assertRaises(FileNotFoundError):
            extraer_vertidas(IFC_BASE.parent / "no_existe.ifc")


class VoxelizacionTests(unittest.TestCase):
    def test_cubo_alineado_reproduce_volumen_y_id(self):
        grilla = voxelizar_vertidas([vertida_cubo()], 0.25)
        self.assertEqual(grilla.forma, (4, 4, 4))
        self.assertEqual(np.unique(grilla.vertida_id).tolist(), [1])
        self.assertAlmostEqual(grilla.volumen_voxelizado_m3, 1.0)
        self.assertAlmostEqual(grilla.error_volumen_relativo_total, 0.0)

    def test_voxelizacion_ifc_aproxima_volumen_global(self):
        vertidas = extraer_vertidas(IFC_BASE)
        grilla_gruesa = voxelizar_vertidas(vertidas, 0.5)
        grilla_fina = voxelizar_vertidas(vertidas, 0.25)
        self.assertEqual(set(np.unique(grilla_fina.vertida_id)) - {0}, {1, 2, 3, 4, 5})
        self.assertLess(grilla_gruesa.error_volumen_relativo_total, 0.02)
        self.assertLess(grilla_fina.error_volumen_relativo_total, 0.02)
        self.assertLess(grilla_fina.error_volumen_relativo_total,
                        grilla_gruesa.error_volumen_relativo_total)

    def test_rechaza_paso_no_positivo(self):
        with self.assertRaises(ValueError):
            voxelizar_vertidas([vertida_cubo()], 0.0)

    def test_detecta_superficies_exteriores_e_interfaz_entre_dos_vertidas(self):
        ids = np.array([1, 2], dtype=np.int32).reshape(2, 1, 1)
        caras = caras_superficie_exterior(ids)
        self.assertEqual(sum(int(np.count_nonzero(c)) for c in caras.values()), 10)
        self.assertEqual(area_interfaces_entre_vertidas(ids, 1.0), {(1, 2): 1.0})

    def test_una_vertida_inactiva_se_convierte_en_vacio_para_la_mascara(self):
        ids = np.array([1, 2], dtype=np.int32).reshape(2, 1, 1)
        activos = np.array([True, False]).reshape(2, 1, 1)
        caras = caras_superficie_exterior(ids, activos)
        self.assertEqual(sum(int(np.count_nonzero(c)) for c in caras.values()), 6)
        self.assertEqual(area_interfaces_entre_vertidas(ids, 1.0, activos), {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
