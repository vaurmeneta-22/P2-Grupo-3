"""Pruebas E11: exportación VTK, snapshot con inicio automático y descarga histórica."""
from datetime import datetime, timedelta, timezone
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.visualization.vtk_export import (
    CODIGO_CARA_RECONSTRUIDA, CODIGO_NUCLEO, Marca, escribir_estado, escribir_marcas, exportar_estados)
from mcplan.weather.loader import construir_aire
from mcplan.weather.openmeteo import construir_url_historica, descargar_historica_y_guardar
from mcplan.weather.serie import SerieTemperatura, a_utc
from mcplan.weather.sitio import Sitio, SitioSinCoordenadas
from mcplan.weather.snapshots import guardar_snapshot, leer_snapshot

H = 3600.0


def grilla_prueba():
    ids = np.zeros((3, 2, 2), dtype=np.int32)
    ids[:, :, 0] = 1
    ids[:2, :, 1] = 2
    return GrillaVoxel((10.0, 20.0, 0.0), (0.5, 0.5, 0.25), ids, {1: 1.5, 2: 0.5})


def leer_vtu(ruta):
    raiz = ET.parse(ruta).getroot()
    pieza = raiz.find("UnstructuredGrid/Piece")

    def arreglo(padre, nombre):
        for nodo in pieza.findall(f"{padre}/DataArray"):
            if nodo.get("Name") == nombre:
                return np.array(nodo.text.split(), dtype=float)
        raise KeyError(nombre)

    campo = {n.get("Name"): float(n.text.split()[0])
             for n in raiz.findall("UnstructuredGrid/FieldData/DataArray")}
    return pieza, arreglo, campo


class ExportacionVtkTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self.grilla = grilla_prueba()
        self.T = np.arange(12, dtype=float).reshape(3, 2, 2) + 20.0
        self.activa = self.grilla.vertida_id == 1          # solo la capa inferior está colocada

    def test_solo_se_exportan_voxeles_activos_con_unidades_y_tiempo(self):
        n = escribir_estado(self.dir / "e.vtu", self.grilla, self.T, self.activa, 36 * H)
        self.assertEqual(n, 6)
        pieza, arreglo, campo = leer_vtu(self.dir / "e.vtu")
        self.assertEqual(int(pieza.get("NumberOfCells")), 6)
        self.assertTrue(np.all(arreglo("Cells", "types") == 12))                 # hexaedros
        self.assertTrue(np.all(arreglo("CellData", "mascara_activa") == 1))
        self.assertTrue(np.all(arreglo("CellData", "vertida_id") == 1))           # la vertida 2 aún no existe
        self.assertEqual(sorted(arreglo("CellData", "temperatura_C")), sorted(self.T[:, :, 0].ravel()))
        self.assertEqual(campo["TimeValue"], 36.0)
        self.assertEqual(campo["tiempo_s"], 36 * H)
        # Coordenadas en metros: origen (10, 20, 0), voxel 0.5 x 0.5 x 0.25, tres voxeles en x.
        coord = arreglo("Points", "coordenadas_m").reshape(-1, 3)
        self.assertAlmostEqual(coord[:, 0].min(), 10.0)
        self.assertAlmostEqual(coord[:, 0].max(), 11.5)
        self.assertAlmostEqual(coord[:, 2].max(), 0.25)

    def test_los_datos_de_campo_declaran_su_numero_de_tuplas(self):
        """Sin NumberOfTuples VTK lee TimeValue/tiempo_h vacíos: hay que declararlo."""
        escribir_estado(self.dir / "t.vtu", self.grilla, self.T, self.activa, 5 * H)
        raiz = ET.parse(self.dir / "t.vtu").getroot()
        arreglos = raiz.findall("UnstructuredGrid/FieldData/DataArray")
        self.assertEqual({a.get("Name") for a in arreglos}, {"TimeValue", "tiempo_h", "tiempo_s"})
        self.assertTrue(all(a.get("NumberOfTuples") == "1" for a in arreglos))

    @unittest.skipUnless(importlib.util.find_spec("pyvista"), "PyVista no instalado")
    def test_pyvista_lee_celdas_marcas_y_tiempo(self):
        import pyvista as pv
        marcas = [Marca(2, CODIGO_NUCLEO, (10.25, 20.25, 0.125), 38.0)]
        escribir_estado(self.dir / "e.vtu", self.grilla, self.T, self.activa, 36 * H)
        escribir_marcas(self.dir / "m.vtu", marcas, 36 * H)
        malla = pv.read(str(self.dir / "e.vtu"))
        self.assertEqual(malla.n_cells, 6)
        self.assertEqual(float(malla.field_data["tiempo_h"][0]), 36.0)
        self.assertEqual(float(pv.read(str(self.dir / "m.vtu"))["temperatura_C"][0]), 38.0)

    def test_caras_expuestas_y_modo_con_voxeles_inactivos(self):
        escribir_estado(self.dir / "a.vtu", self.grilla, self.T, self.activa, 0.0)
        _, arreglo, _ = leer_vtu(self.dir / "a.vtu")
        self.assertTrue(np.all(arreglo("CellData", "n_caras_expuestas") >= 1))
        n = escribir_estado(self.dir / "b.vtu", self.grilla, self.T, self.activa, 0.0, solo_activas=False)
        self.assertEqual(n, 10)
        _, arreglo, _ = leer_vtu(self.dir / "b.vtu")
        self.assertEqual(int(arreglo("CellData", "mascara_activa").sum()), 6)

    def test_marcas_distinguen_voxel_y_cara_reconstruida(self):
        marcas = [Marca(2, CODIGO_NUCLEO, (10.25, 20.25, 0.125), 38.0),
                  Marca(2, CODIGO_CARA_RECONSTRUIDA, (11.5, 20.25, 0.125), 24.0)]
        escribir_marcas(self.dir / "m.vtu", marcas, 12 * H)
        pieza, arreglo, campo = leer_vtu(self.dir / "m.vtu")
        self.assertEqual(int(pieza.get("NumberOfPoints")), 2)
        self.assertTrue(np.all(arreglo("Cells", "types") == 1))
        self.assertEqual(list(arreglo("PointData", "codigo_marca")), [1, 3])
        self.assertEqual(list(arreglo("PointData", "temperatura_C")), [38.0, 24.0])
        self.assertEqual(campo["tiempo_h"], 12.0)
        with self.assertRaises(ValueError):
            escribir_marcas(self.dir / "x.vtu", [], 0.0)

    def test_pvd_y_estados_mas_cercanos(self):
        casos = self

        class Resultado:
            tiempos_s = np.array([0.0, 12 * H, 24 * H, 48 * H])
            temperaturas_c = tuple(casos.T + i for i in range(4))
            mascaras_activas = tuple(casos.activa for _ in range(4))

        resumen = exportar_estados(self.dir / "vtk", self.grilla, Resultado, [0.0, 13.0, 50.0],
                                   marcas_por_estado=lambda n, t: [Marca(1, 1, (10.2, 20.2, 0.1), 30.0)])
        self.assertEqual([round(e["tiempo_h"]) for e in resumen["estados"]], [0, 12, 48])
        raiz = ET.parse(self.dir / "vtk" / "estado.pvd").getroot()
        self.assertEqual([float(d.get("timestep")) for d in raiz.iter("DataSet")], [0.0, 12.0, 48.0])
        for d in raiz.iter("DataSet"):
            self.assertTrue((self.dir / "vtk" / d.get("file")).is_file())
        self.assertTrue((self.dir / "vtk" / "marcas.pvd").is_file())
        leyenda = json.loads((self.dir / "vtk" / "leyenda_marcas.json").read_text(encoding="utf-8"))
        self.assertIn("3", leyenda["codigo_marca"])


class SnapshotE11Test(unittest.TestCase):
    RESPUESTA = json.dumps({"utc_offset_seconds": 0, "hourly_units": {"temperature_2m": "°C"},
                            "hourly": {"time": ["2025-11-11T00:00", "2025-11-11T01:00"],
                                       "temperature_2m": [14.0, 13.0]}}).encode("utf-8")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.raiz = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def snapshot(self, horas=72, inicio="2026-11-11T12:00:00Z"):
        t0 = a_utc(inicio)
        tiempos = tuple(t0 + timedelta(hours=h) for h in range(horas + 1))
        serie = SerieTemperatura(tiempos, np.linspace(10.0, 20.0, horas + 1),
                                 tuple("pronostico" for _ in tiempos))
        return guardar_snapshot(self.raiz / "data" / "weather", bruto=b"{}",
                                nombre_original="original.json", serie=serie, fuente="prueba",
                                modelo="m", latitud=-30.6, longitud=-71.2,
                                descargado_utc=datetime(2026, 11, 10, tzinfo=timezone.utc))

    def config(self, carpeta):
        return {"fuente": "snapshot", "snapshot": {"carpeta": str(carpeta.relative_to(self.raiz))}}

    def test_inicio_automatico_usa_el_primer_instante_del_snapshot(self):
        carpeta = self.snapshot()
        aire = construir_aire(self.config(carpeta), tiempo_final_s=48 * H, raiz_repo=self.raiz)
        self.assertEqual(aire.descripcion["inicio_obra_utc"], "2026-11-11T12:00:00Z")
        self.assertAlmostEqual(aire.funcion(0.0), 10.0)
        self.assertAlmostEqual(aire.funcion(72 * H), 20.0)
        self.assertTrue(aire.cobertura.cubre_completo)
        with self.assertRaises(ValueError) as error:       # horizonte mayor que el intervalo cubierto
            construir_aire(self.config(carpeta), tiempo_final_s=96 * H, raiz_repo=self.raiz)
        self.assertIn("2026-11-14T12:00:00Z", str(error.exception))

    def test_el_loader_informa_tipo_y_metadatos_del_snapshot(self):
        carpeta = self.snapshot()
        aire = construir_aire(self.config(carpeta), tiempo_final_s=H, raiz_repo=self.raiz)
        self.assertFalse(aire.es_sintetica)
        self.assertEqual(aire.descripcion["tipo"], "pronostico")
        self.assertEqual(leer_snapshot(carpeta).metadatos["latitude"], -30.6)

    def test_descarga_historica_es_reanalisis_y_se_bloquea_sin_coordenadas(self):
        url = construir_url_historica(-30.6, -71.2, "2025-11-11", "2025-11-25")
        self.assertIn("archive-api.open-meteo.com", url)
        self.assertIn("start_date=2025-11-11", url)
        sitio = Sitio(3, "Ovalle / Ecos del Mar", "Coquimbo", -30.6, -71.2, "prueba")
        carpeta = descargar_historica_y_guardar(
            sitio, self.raiz / "w", "2025-11-11", "2025-11-25", obtener=lambda _u: self.RESPUESTA,
            ahora_utc=datetime(2026, 11, 10, tzinfo=timezone.utc))
        meta = json.loads((carpeta / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta["series_kind"], "reanalisis")
        self.assertEqual(meta["coordinates_source"], "prueba")          # procedencia de las coordenadas
        self.assertEqual((meta["latitude"], meta["longitude"]), (-30.6, -71.2))
        self.assertIsNone(meta["issued_at"])
        llamadas = []
        sin = Sitio(3, "Ovalle / Ecos del Mar", "Coquimbo", None, None, "pendiente de sitios.yaml")
        with self.assertRaises(SitioSinCoordenadas):
            descargar_historica_y_guardar(sin, self.raiz, "2025-11-11", "2025-11-25",
                                          obtener=lambda u: llamadas.append(u) or b"")
        self.assertEqual(llamadas, [])


if __name__ == "__main__":
    unittest.main()
