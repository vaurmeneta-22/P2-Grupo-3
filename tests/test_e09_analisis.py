"""Pruebas de las funciones de análisis y figuras E09."""
from __future__ import annotations

import csv
import dataclasses
import sys
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

from mcplan.geometry.ifc import VertidaIFC, extraer_vertidas
from mcplan.geometry.voxelizacion import (
    GrillaVoxel,
    caras_superficie_exterior,
)
from mcplan.thermal.analisis import (
    corte_termico,
    cortes_por_punto,
    eventos_activacion,
    metricas_tmax,
    seleccionar_puntos_representativos,
    serie_punto,
)
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.vertidas import ResultadoEtapas
from scripts.demo_pipeline_ifc_3d import (
    COLORES_VERTIDA_VISTAS,
    _offsets_explotados,
    calcular_perfil_vertical_nucleo_superficie,
    dibujar_cortes_termicos,
    dibujar_geometria_ifc,
    dibujar_geometria_ifc_ensamblada,
    dibujar_geometria_ifc_explotada,
    dibujar_historial,
    dibujar_puntos_representativos,
    guardar_metricas_e09,
    guardar_perfil_vertical_e09,
    guardar_puntos_e09,
    triangulos_exteriores_ensamblados,
)


def crear_grilla() -> GrillaVoxel:
    ids = np.zeros((3, 3, 4), dtype=np.int32)
    ids[:, :, :2] = 1
    ids[:, :, 2:] = 2
    ids[2, 2, :] = 0
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 3.0, 2: 3.0})


def crear_resultado(grilla: GrillaVoxel) -> ResultadoEtapas:
    ids = grilla.vertida_id
    activa_1 = ids == 1
    activa_final = ids > 0

    campo_0 = np.zeros(ids.shape, dtype=float)
    campo_0[activa_1] = 20.0

    campo_1 = np.zeros(ids.shape, dtype=float)
    campo_1[ids == 1] = 20.0
    campo_1[ids == 2] = 25.0
    campo_1[1, 1, 1] = 50.0   # P1: interior de la primera vertida.
    campo_1[0, 1, 2] = 45.0   # P2: superficie expuesta más caliente.
    campo_1[1, 1, 2] = 40.0   # P3: interior de la última vertida.

    campo_2 = np.zeros(ids.shape, dtype=float)
    campo_2[activa_final] = 30.0
    campo_2[0, 1, 2] = 32.0

    tiempos = np.array([0.0, 6.0 * 3600.0, 12.0 * 3600.0])
    campos = (campo_0, campo_1, campo_2)
    mascaras = (activa_1.copy(), activa_final.copy(), activa_final.copy())
    maximas = np.array([np.max(c[m]) for c, m in zip(campos, mascaras)])
    promedios = np.array([np.mean(c[m]) for c, m in zip(campos, mascaras)])
    return ResultadoEtapas(
        tiempos_s=tiempos,
        temperaturas_c=campos,
        mascaras_activas=mascaras,
        ids_activados=((1,), (2,), ()),
        temperaturas_maximas_c=maximas,
        temperaturas_promedio_c=promedios,
        energia_fuente_j=0.0,
        energia_activacion_j=0.0,
        energia_frontera_j=0.0,
        energia_final_j=0.0,
        error_balance_j=0.0,
        ultimo_paso_maximo_estable_s=600.0,
    )


IFC_BASE = RAIZ / "models" / "ifc" / "Fundacion_1_con_propiedades.ifc"


def vertida_cubo(lado_m: float = 1.0) -> VertidaIFC:
    vertices = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float) * lado_m
    triangulos = np.array([
        [0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
        [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
        [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7],
    ], dtype=np.int64)
    return VertidaIFC(
        1, "Vertida_01", "cubo", 0.0, 18.0, "M1",
        vertices, triangulos,
        {"VertidaID": 1, "Inicio_h": 0.0,
         "TemperaturaColocacion_C": 18.0, "TipoHormigon": "M1"},
    )


class AnalisisE09Tests(unittest.TestCase):
    def setUp(self):
        self.grilla = crear_grilla()
        self.resultado = crear_resultado(self.grilla)
        self.puntos = seleccionar_puntos_representativos(self.grilla, self.resultado)

    def test_serie_punto_respeta_activacion_y_valores(self):
        horas, temperaturas = serie_punto(self.resultado, (1, 1, 2))
        np.testing.assert_allclose(horas, [0.0, 6.0, 12.0])
        self.assertTrue(np.isnan(temperaturas[0]))
        self.assertAlmostEqual(temperaturas[1], 40.0)
        self.assertAlmostEqual(temperaturas[2], 30.0)

        _, temperaturas_p1 = serie_punto(self.resultado, (1, 1, 1))
        np.testing.assert_allclose(temperaturas_p1, [20.0, 50.0, 30.0])

    def test_serie_punto_rechaza_ijk_fuera_de_grilla(self):
        with self.assertRaises(ValueError):
            serie_punto(self.resultado, (3, 0, 0))

    def test_selecciona_tres_puntos_con_categorias_esperadas(self):
        self.assertEqual([p.etiqueta for p in self.puntos], ["P1", "P2", "P3"])
        self.assertEqual([p.categoria for p in self.puntos],
                         ["pico_global", "superficie_expuesta",
                          "nucleo_ultima_vertida"])
        self.assertEqual([p.ijk for p in self.puntos],
                         [(1, 1, 1), (0, 1, 2), (1, 1, 2)])
        self.assertEqual([p.vertida_id for p in self.puntos], [1, 2, 2])
        np.testing.assert_allclose(self.puntos[0].xyz_m, (1.5, 1.5, 1.5))

    def test_p2_esta_en_superficie_y_p3_es_interior_de_la_ultima_vertida(self):
        activa = self.resultado.mascaras_activas[-1]
        caras = caras_superficie_exterior(self.grilla.vertida_id, activa)
        superficie = np.logical_or.reduce(tuple(caras.values()))
        ijk_p2 = self.puntos[1].ijk
        ijk_p3 = self.puntos[2].ijk
        self.assertTrue(superficie[ijk_p2])
        self.assertFalse(superficie[ijk_p3])
        self.assertEqual(int(self.grilla.vertida_id[ijk_p3]), 2)

    def test_seleccion_es_determinista(self):
        repetidos = seleccionar_puntos_representativos(self.grilla, self.resultado)
        self.assertEqual(self.puntos, repetidos)

    def test_corte_termico_usa_el_tiempo_mas_cercano_y_enmascara_exterior(self):
        corte = corte_termico(self.grilla, self.resultado, 5.2 * 3600.0, "z", 2)
        self.assertEqual(corte.eje, "z")
        self.assertEqual(corte.indice, 2)
        self.assertAlmostEqual(corte.tiempo_h, 6.0)
        self.assertAlmostEqual(corte.coordenada_m, 2.5)
        self.assertEqual(corte.temperaturas_c.shape, (3, 3))
        self.assertAlmostEqual(corte.temperaturas_c[1, 1], 40.0)
        self.assertAlmostEqual(corte.temperaturas_c[0, 1], 45.0)
        self.assertTrue(np.isnan(corte.temperaturas_c[2, 2]))
        np.testing.assert_allclose(corte.extension_m, (3.0, 3.0))

    def test_tres_cortes_pasan_por_el_mismo_punto(self):
        cortes = cortes_por_punto(self.grilla, self.resultado, (1, 1, 2), 6 * 3600.0)
        self.assertEqual([c.eje for c in cortes], ["x", "y", "z"])
        self.assertEqual([c.indice for c in cortes], [1, 1, 2])
        self.assertEqual(cortes[0].temperaturas_c.shape, (3, 4))
        self.assertEqual(cortes[1].temperaturas_c.shape, (3, 4))

    def test_corte_rechaza_eje_invalido(self):
        with self.assertRaises(ValueError):
            corte_termico(self.grilla, self.resultado, 0.0, "w", 1)

    def test_metricas_tmax_y_eventos(self):
        met = metricas_tmax(self.resultado)
        self.assertAlmostEqual(met.t_pico_h, 6.0)
        self.assertAlmostEqual(met.tmax_pico_c, 50.0)
        self.assertAlmostEqual(met.tmax_final_c, 32.0)
        self.assertEqual(met.n_muestras, 3)
        self.assertEqual(eventos_activacion(self.resultado),
                         ((0.0, (1,)), (6.0, (2,))))

    def test_perfil_vertical_separa_nucleo_y_superficie(self):
        campos = [campo.copy() for campo in self.resultado.temperaturas_c]
        campos[1][0, 0, 2] = 30.0   # núcleo interior en la columna (0, 0)
        campos[1][0, 0, 3] = 22.0   # superficie superior en la misma columna
        resultado = dataclasses.replace(
            self.resultado, temperaturas_c=tuple(campos))
        filas, columna = calcular_perfil_vertical_nucleo_superficie(
            self.grilla, resultado, ijk_referencia=(0, 0))
        self.assertEqual(columna, (0, 0))
        self.assertEqual(len(filas), 2)  # t = 6 h y t = 12 h con 4 capas activas
        self.assertAlmostEqual(filas[0][4], 30.0)
        self.assertAlmostEqual(filas[0][6], 22.0)
        self.assertAlmostEqual(filas[0][7], 8.0)
        for fila in filas:
            self.assertAlmostEqual(fila[7], fila[4] - fila[6])

    def test_exporta_perfil_vertical_con_delta(self):
        with tempfile.TemporaryDirectory() as temporal:
            ruta = Path(temporal) / "perfil.csv"
            filas = [
                (6.0, 0.0, 0.0, 2.0, 30.0, 3.0, 22.0, 8.0),
                (12.0, 0.0, 0.0, 2.0, 31.0, 3.0, 25.0, 6.0),
            ]
            guardar_perfil_vertical_e09(ruta, filas)
            with ruta.open(encoding="utf-8") as archivo:
                leidas = list(csv.DictReader(archivo))
        self.assertEqual(leidas[0]["delta_c"], "8.00000000")
        self.assertEqual(leidas[1]["t_superficie_c"], "25.00000000")
        self.assertEqual(len(leidas[0]), 8)


class MallasIFCTests(unittest.TestCase):
    def test_cinco_mallas_cerradas_y_volumen_coincide(self):
        vertidas = extraer_vertidas(IFC_BASE)
        self.assertEqual(len(vertidas), 5)
        total_triangulos = 0
        z_anterior = 0.0
        for vertida in vertidas:
            unicos, inverso = np.unique(
                np.round(vertida.vertices_m, 6), axis=0, return_inverse=True
            )
            self.assertEqual(len(unicos), len(vertida.vertices_m))
            aristas = {}
            for triangulo in vertida.triangulos:
                ids = inverso[triangulo].tolist()
                for a, b in ((ids[0], ids[1]), (ids[1], ids[2]),
                             (ids[2], ids[0])):
                    clave = tuple(sorted((a, b)))
                    aristas[clave] = aristas.get(clave, 0) + 1
            self.assertTrue(all(cantidad == 2 for cantidad in aristas.values()))

            a = vertida.vertices_m[vertida.triangulos[:, 0]]
            b = vertida.vertices_m[vertida.triangulos[:, 1]]
            c = vertida.vertices_m[vertida.triangulos[:, 2]]
            volumen = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
            self.assertAlmostEqual(volumen, vertida.volumen_m3, places=6)
            self.assertAlmostEqual(float(vertida.minimo_m[2]), z_anterior, places=6)
            z_anterior = float(vertida.maximo_m[2])
            total_triangulos += len(vertida.triangulos)

        self.assertEqual(total_triangulos, 300)
        self.assertAlmostEqual(z_anterior, 4.3, places=6)


class FigurasE09Tests(unittest.TestCase):
    def setUp(self):
        self.grilla = crear_grilla()
        self.resultado = crear_resultado(self.grilla)
        self.puntos = seleccionar_puntos_representativos(self.grilla, self.resultado)
        self.temporal = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporal.cleanup)
        self.carpeta = Path(self.temporal.name)

    def _comprobar_png(self, figura, ruta: Path):
        plt.close(figura)
        self.assertTrue(ruta.is_file())
        self.assertGreater(ruta.stat().st_size, 1000)

    def _conteo_caras_por_vertida(self, figura):
        colecciones = [c for c in figura.axes[0].collections
                       if hasattr(c, "get_facecolors")]
        self.assertEqual(len(colecciones), 1)
        colores = np.asarray(colecciones[0].get_facecolors())[:, :3]
        bases = np.array([
            matplotlib.colors.to_rgb(COLORES_VERTIDA_VISTAS[vid])
            for vid in range(1, 6)
        ])
        bases = bases / np.linalg.norm(bases, axis=1, keepdims=True)
        normales = colores / np.linalg.norm(colores, axis=1, keepdims=True)
        asignado = (normales @ bases.T).argmax(axis=1)
        conteo = np.bincount(asignado, minlength=5).tolist()
        return conteo, len(colores)

    def test_figura_bim_reutiliza_mallas_del_ifc(self):
        ruta = self.carpeta / "bim.png"
        figura = dibujar_geometria_ifc([vertida_cubo()], ruta)
        self._comprobar_png(figura, ruta)

    def test_figura_bim_con_las_cinco_vertidas_reales(self):
        vertidas = extraer_vertidas(IFC_BASE)
        self.assertEqual(len(vertidas), 5)
        self.assertTrue(all(len(v.triangulos) > 0 for v in vertidas))
        ruta = self.carpeta / "bim_ifc_real.png"
        figura = dibujar_geometria_ifc(vertidas, ruta)
        leyenda = figura.axes[0].get_legend()
        textos = [texto.get_text() for texto in leyenda.get_texts()]
        self.assertEqual(len(textos), 5)
        self.assertTrue(any("inicio 0 h" in texto for texto in textos))
        self.assertTrue(any("inicio 24 h" in texto for texto in textos))
        self.assertTrue(all("triángulos IFC" in texto for texto in textos))
        self._comprobar_png(figura, ruta)

    def test_vista_ensamblada_conserva_envolvente_exterior(self):
        vertidas = extraer_vertidas(IFC_BASE)
        visibles = triangulos_exteriores_ensamblados(vertidas)
        conteo = {vid: len(indices) for vid, indices in visibles.items()}
        self.assertEqual(conteo, {1: 22, 2: 16, 3: 16, 4: 22, 5: 142})
        self.assertEqual(sum(conteo.values()), 218)
        self.assertEqual(sum(len(v.triangulos) for v in vertidas), 300)

        # Las tapas/base internas se omiten, pero ninguna cara lateral.
        for vertida in vertidas:
            triangulos_visibles = visibles[vertida.vertida_id]
            for indice, triangulo in enumerate(vertida.triangulos):
                a, b, c = vertida.vertices_m[triangulo]
                normal = np.cross(b - a, c - a)
                norma = float(np.linalg.norm(normal))
                horizontal = abs(float(normal[2])) >= 0.999 * norma
                if not horizontal:
                    self.assertIn(indice, triangulos_visibles)

        # Ningún plano horizontal se conserva en dos vertidas distintas.
        z_horizontales = {}
        for vertida in vertidas:
            for indice in visibles[vertida.vertida_id]:
                a, b, c = vertida.vertices_m[vertida.triangulos[indice]]
                normal = np.cross(b - a, c - a)
                norma = float(np.linalg.norm(normal))
                if abs(float(normal[2])) < 0.999 * norma:
                    continue
                z = round(float((a[2] + b[2] + c[2]) / 3.0), 6)
                z_horizontales.setdefault(z, set()).add(vertida.vertida_id)
        compartidos = {z: ids for z, ids in z_horizontales.items() if len(ids) > 1}
        self.assertEqual(compartidos, {})

    def test_vista_ensamblada_etiquetas_y_png(self):
        vertidas = extraer_vertidas(IFC_BASE)
        ruta = self.carpeta / "bim_ensamblada.png"
        figura = dibujar_geometria_ifc_ensamblada(vertidas, ruta)
        eje = figura.axes[0]
        conteo, total = self._conteo_caras_por_vertida(figura)
        self.assertEqual(total, 218)
        self.assertEqual(conteo, [22, 16, 16, 22, 142])
        rotulos = [texto.get_text() for texto in eje.texts
                   if texto.get_text().startswith("V")]
        self.assertEqual(len(rotulos), 5)
        for vid, hora in ((1, "0 h"), (2, "6 h"), (3, "12 h"),
                          (4, "18 h"), (5, "24 h")):
            self.assertTrue(any(f"V{vid} · inicio {hora}" in rotulo
                                for rotulo in rotulos))
        leyenda = [texto.get_text() for texto in eje.get_legend().get_texts()]
        self.assertEqual(len(leyenda), 5)
        self.assertTrue(all("inicio" in texto for texto in leyenda))
        self.assertIn("218 de 300", figura.texts[0].get_text())
        self._comprobar_png(figura, ruta)

    def test_vista_explotada_desplaza_solo_el_dibujo(self):
        vertidas = extraer_vertidas(IFC_BASE)
        originales = [(v.vertices_m.copy(), v.volumen_m3,
                       v.minimo_m.copy(), v.maximo_m.copy())
                      for v in vertidas]
        offsets = _offsets_explotados(vertidas, 0.85)
        self.assertEqual(offsets, {1: 0.0, 2: 0.85, 3: 1.70,
                                   4: 2.55, 5: 3.40})

        ruta = self.carpeta / "bim_explotada.png"
        figura = dibujar_geometria_ifc_explotada(vertidas, ruta, 0.85)
        eje = figura.axes[0]
        conteo, total = self._conteo_caras_por_vertida(figura)
        self.assertEqual(total, 300)
        self.assertEqual(conteo, [28, 28, 28, 28, 188])
        rotulos = [texto.get_text() for texto in eje.texts
                   if texto.get_text().startswith("V")]
        self.assertEqual(len(rotulos), 5)
        for vid, hora in ((1, "0 h"), (2, "6 h"), (3, "12 h"),
                          (4, "18 h"), (5, "24 h")):
            self.assertTrue(any(f"V{vid} · inicio {hora}" in rotulo
                                for rotulo in rotulos))
        self.assertTrue(any("explotada" in texto.get_text()
                            for texto in [eje.title]))
        self.assertTrue(any("Visualización explicativa" in texto.get_text()
                            for texto in figura.texts))
        self.assertTrue(any("no geometría de simulación" in texto.get_text()
                            for texto in figura.texts))

        # Cada vertida conserva vértices, volumen y cajas originales.
        for vertida, (vertices, volumen, minimo, maximo) in zip(vertidas, originales):
            np.testing.assert_allclose(vertida.vertices_m, vertices)
            self.assertAlmostEqual(vertida.volumen_m3, volumen, places=9)
            np.testing.assert_allclose(vertida.minimo_m, minimo)
            np.testing.assert_allclose(vertida.maximo_m, maximo)
        self._comprobar_png(figura, ruta)

    def test_figura_cortes_incluye_titulos_y_contexto(self):
        ruta = self.carpeta / "cortes_contexto.png"
        figura = dibujar_cortes_termicos(
            self.grilla, self.resultado, (1, 1, 2),
            [6 * 3600.0, 12 * 3600.0], ruta,
        )
        titulos = [eje.get_title() for eje in figura.axes]
        self.assertTrue(any("Corte vertical YZ" in titulo for titulo in titulos))
        self.assertTrue(any("Corte vertical XZ" in titulo for titulo in titulos))
        self.assertTrue(any("Corte horizontal XY" in titulo for titulo in titulos))
        textos = [texto.get_text() for texto in figura.texts]
        self.assertTrue(any("máximo térmico global" in texto for texto in textos))
        self.assertTrue(any("final de la simulación" in texto for texto in textos))
        self.assertTrue(any("exteriores al hormigón" in texto for texto in textos))
        self._comprobar_png(figura, ruta)

    def test_figura_tres_puntos_separa_maximos_y_criterios(self):
        ruta = self.carpeta / "puntos_contexto.png"
        figura = dibujar_puntos_representativos(
            self.grilla, self.resultado, self.puntos, ruta,
        )
        panel = figura.axes[1]
        self.assertTrue(panel.tables)
        tabla = panel.tables[0]
        filas = [
            [tabla[(fila, columna)].get_text().get_text() for columna in range(4)]
            for fila in range(4)
        ]
        self.assertEqual([fila[0] for fila in filas[1:]], ["P1", "P2", "P3"])
        self.assertEqual(filas[1][1], "50.00")
        self.assertTrue(any("Criterio de selección" in texto.get_text()
                            for texto in panel.texts))
        self.assertTrue(any("Punto más caliente" in texto.get_text()
                            for texto in panel.texts))
        self._comprobar_png(figura, ruta)

    def test_figura_cortes_termicos(self):
        ruta = self.carpeta / "cortes.png"
        figura = dibujar_cortes_termicos(
            self.grilla, self.resultado, (1, 1, 2),
            [6 * 3600.0, 12 * 3600.0], ruta,
        )
        self._comprobar_png(figura, ruta)

    def test_figura_tres_puntos(self):
        ruta = self.carpeta / "puntos.png"
        figura = dibujar_puntos_representativos(
            self.grilla, self.resultado, self.puntos, ruta,
        )
        self._comprobar_png(figura, ruta)

    def test_figura_tmax_con_pico(self):
        ruta = self.carpeta / "tmax.png"
        figura = dibujar_historial(self.resultado, ruta)
        self._comprobar_png(figura, ruta)

    def test_exporta_metricas_y_puntos(self):
        ruta_metricas = self.carpeta / "metricas.csv"
        ruta_puntos = self.carpeta / "puntos.csv"
        ruta_series = self.carpeta / "series.csv"
        hidratacion = ModeloPulsoCosenoidal(500.0, 48 * 3600.0)
        guardar_metricas_e09(
            ruta_metricas, grilla=self.grilla, resultado=self.resultado,
            puntos=self.puntos, hidratacion=hidratacion,
            config={"condicion_superficie": {
                "tipo": "robin", "coeficiente_h_w_m2_k": 8.0,
                "temperatura_aire_c": 20.0,
            }},
            dx=1.0, area_exterior_m2=12.0, interfaces={(1, 2): 9.0},
        )
        guardar_puntos_e09(ruta_puntos, self.resultado, self.puntos,
                           ruta_series=ruta_series)
        with ruta_metricas.open(encoding="utf-8") as archivo:
            filas = list(csv.DictReader(archivo))
        valores = {fila["metrica"]: fila["valor"] for fila in filas}
        self.assertEqual(valores["tmax_pico_h"], "6.000000")
        self.assertEqual(valores["condicion_superficie"], "robin")
        with ruta_puntos.open(encoding="utf-8") as archivo:
            filas = list(csv.DictReader(archivo))
        self.assertEqual([fila["etiqueta"] for fila in filas], ["P1", "P2", "P3"])
        with ruta_series.open(encoding="utf-8") as archivo:
            series = list(csv.DictReader(archivo))
        self.assertEqual(len(series), 9)
        p1 = [fila for fila in series if fila["etiqueta"] == "P1"]
        self.assertEqual([fila["temperatura_c"] for fila in p1],
                         ["20.00000000", "50.00000000", "30.00000000"])
        p3 = [fila for fila in series if fila["etiqueta"] == "P3"]
        self.assertEqual(p3[0]["temperatura_c"], "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
