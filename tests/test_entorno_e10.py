"""Pruebas E10: suelo efectivo, resistencias en serie, retiro y cobertura del loader."""
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
import tempfile
import unittest

import numpy as np

from mcplan.geometry.ifc import VertidaIFC
from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.bordes import (
    ModeloBordes, ModeloSuelo, Proteccion, crear_modelo_bordes, proteccion_desde_config,
    resistencia_de_capa)
from mcplan.thermal.diffusion_nd import DifusionCartesiana, robin
from mcplan.thermal.extremos import SeguimientoExtremos
from mcplan.thermal.hidratacion import ModeloPulsoCosenoidal
from mcplan.thermal.seguimiento import temperatura_cara_estimada
from mcplan.thermal.vertidas import simular_vertidas
from mcplan.weather.loader import construir_aire
from mcplan.weather.openmeteo import construir_url, descargar_y_guardar, parsear_respuesta
from mcplan.weather.serie import SerieTemperatura, a_utc, normalizar_serie
from mcplan.weather.sintetica import serie_sintetica
from mcplan.weather.sitio import Sitio, SitioSinCoordenadas, sitio_desde_config
from mcplan.weather.snapshots import (
    cobertura, extender_serie, guardar_snapshot, leer_snapshot, listar_snapshots)

H = 3600.0
CARAS = ("x-", "x+", "y-", "y+", "z-", "z+")


def vertida(identificador, inicio_h=0.0, temperatura_c=30.0, propiedades=None):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    triangulos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    return VertidaIFC(identificador, f"Vertida_{identificador:02d}", "gid", inicio_h,
                      temperatura_c, "M1", vertices, triangulos, propiedades or {})


def grilla_bloque(n=2):
    """Un bloque n x n x n de 1 m, una sola vertida (ID 1)."""
    ids = np.ones((n, n, n), dtype=np.int32)
    return GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: float(n ** 3)})


def solver_bloque(grilla, **extra):
    mascara = grilla.vertida_id > 0
    return DifusionCartesiana(
        forma=grilla.forma, dx_m=grilla.paso_m, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
        calor_especifico_j_kg_k=900.0, fronteras={c: robin(10.0, 20.0) for c in CARAS},
        mascara_activa=mascara, **extra)


class SueloTest(unittest.TestCase):
    def test_suelo_efectivo_es_robin_con_h_eq_k_sobre_L(self):
        suelo = ModeloSuelo(temperatura_c=10.0, conductividad_w_m_k=1.5, espesor_efectivo_m=0.75)
        self.assertAlmostEqual(suelo.h_equivalente_w_m2_k, 2.0)
        self.assertEqual(suelo.condicion().tipo, "robin")
        with self.assertRaises(ValueError):
            ModeloSuelo(10.0, 0.0, 1.0)

    def test_flujo_de_la_base_sigue_la_formula_en_serie_y_las_caras_laterales_usan_aire(self):
        grilla = grilla_bloque()
        suelo = ModeloSuelo(10.0, 1.5, 0.75)             # R_suelo = 0.5 m² K/W
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), suelo=suelo)
        solver = solver_bloque(grilla, fronteras_celda=bordes.fronteras_celda())
        T = np.full(grilla.forma, 30.0)
        flujos = solver.flujos_frontera_w(T, 0.0)
        area, d, k = 1.0, 1.0, 2.0
        esperado_base = area * (10.0 - 30.0) / (d / (2 * k) + 0.5)
        esperado_aire = area * (20.0 - 30.0) / (d / (2 * k) + 1.0 / 10.0)
        self.assertTrue(np.allclose(flujos["z-"][:, :, 0], esperado_base))
        self.assertTrue(np.allclose(flujos["x+"][-1], esperado_aire))
        self.assertTrue(np.allclose(flujos["z+"][:, :, -1], esperado_aire))

    def test_solo_la_base_toca_el_suelo(self):
        grilla = grilla_bloque(3)
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), suelo=ModeloSuelo(10.0, 1.5, 0.75))
        self.assertEqual(int(bordes.mascara_base.sum()), 9)
        self.assertTrue(np.all(bordes.mascara_base[:, :, 0]))
        self.assertFalse(np.any(bordes.mascara_base[:, :, 1:]))


class BaseGeometricaTest(unittest.TestCase):
    def test_voladizo_con_cara_hacia_abajo_no_toca_el_suelo(self):
        # Corte x-z: base en k=0 (x=0..2) y un voladizo en k=1 que sobresale en x=3 sin apoyo.
        ids = np.zeros((4, 1, 2), dtype=np.int32)
        ids[0:3, 0, 0] = 1
        ids[0:4, 0, 1] = 1
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 7.0})
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), suelo=ModeloSuelo(10.0, 1.5, 0.75))
        self.assertEqual(int(bordes.mascara_base.sum()), 3)
        self.assertTrue(np.all(np.argwhere(bordes.mascara_base)[:, 2] == 0))
        solver = DifusionCartesiana(
            forma=grilla.forma, dx_m=grilla.paso_m, conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
            calor_especifico_j_kg_k=900.0, fronteras={c: robin(10.0, 20.0) for c in CARAS},
            mascara_activa=ids > 0, fronteras_celda=bordes.fronteras_celda())
        flujos = solver.flujos_frontera_w(np.full(grilla.forma, 30.0), 0.0)["z-"]
        self.assertAlmostEqual(flujos[0, 0, 0], (10.0 - 30.0) / (0.25 + 0.5))      # base: suelo
        self.assertAlmostEqual(flujos[3, 0, 1], (20.0 - 30.0) / (0.25 + 0.1))      # voladizo: aire

    def test_la_base_se_ubica_por_la_capa_mas_baja_aunque_no_sea_el_indice_cero(self):
        ids = np.zeros((2, 2, 4), dtype=np.int32)
        ids[:, :, 1:3] = 1
        grilla = GrillaVoxel((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), ids, {1: 8.0})
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), suelo=ModeloSuelo(10.0, 1.5, 0.75))
        self.assertEqual(bordes.indice_base, 1)
        self.assertEqual(int(bordes.mascara_base.sum()), 4)


class ResistenciasTest(unittest.TestCase):
    def test_capas_en_serie_se_suman_con_media_celda_y_conveccion(self):
        grilla = grilla_bloque()
        r_moldaje = resistencia_de_capa({"espesor_m": 0.02, "conductividad_w_m_k": 0.1})    # 0.2
        r_aislacion = resistencia_de_capa({"resistencia_m2_k_w": 0.75})
        protecciones = (Proteccion(1, "moldaje", ("x+",), r_moldaje, None),
                        Proteccion(1, "aislacion", ("x+",), r_aislacion, None))
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), protecciones=protecciones)
        solver = solver_bloque(grilla, resistencias_borde=bordes.resistencias_en(0.0))
        T = np.full(grilla.forma, 30.0)
        flujos = solver.flujos_frontera_w(T, 0.0)
        r_total = 1.0 / (2 * 2.0) + 0.2 + 0.75 + 1.0 / 10.0
        self.assertTrue(np.allclose(flujos["x+"][-1], (20.0 - 30.0) / r_total))
        # Una cara sin protección conserva la resistencia sin capas.
        self.assertTrue(np.allclose(flujos["x-"][0], (20.0 - 30.0) / (0.25 + 0.1)))

    def test_la_cara_estimada_del_hormigon_incluye_la_proteccion_en_serie(self):
        cond = robin(10.0, 20.0)
        sin = temperatura_cara_estimada(cond, 30.0, 2.0, 1.0, 0.0)
        con = temperatura_cara_estimada(cond, 30.0, 2.0, 1.0, 0.0, resistencia_extra_m2_k_w=0.95)
        self.assertAlmostEqual(sin, 30.0 + (20.0 - 30.0) / (0.25 + 0.1) * 0.25)
        self.assertAlmostEqual(con, 30.0 + (20.0 - 30.0) / (0.25 + 0.95 + 0.1) * 0.25)
        self.assertGreater(con, sin)        # con protección la cara del hormigón queda más caliente

    def test_el_solver_y_la_cara_estimada_usan_la_misma_condicion(self):
        """q = (T_aire - T_cara)/(R'' + 1/h) con la T_cara estimada debe igualar el flujo del solver."""
        grilla = grilla_bloque()
        bordes = ModeloBordes(grilla, robin(10.0, 20.0),
                              protecciones=(Proteccion(1, "m", ("x+",), 0.9, None),))
        solver = solver_bloque(grilla, resistencias_borde=bordes.resistencias_en(0.0))
        T = np.full(grilla.forma, 35.0)
        q_solver = solver.flujos_frontera_w(T, 0.0)["x+"][-1, 0, 0]
        cond, r_extra = solver.parametros_borde_celda("x+", (1, 0, 0))
        t_cara = temperatura_cara_estimada(cond, 35.0, 2.0, 1.0, 0.0, r_extra)
        self.assertAlmostEqual(q_solver, (20.0 - t_cara) / (0.9 + 1.0 / 10.0))
        self.assertAlmostEqual(q_solver, (t_cara - 35.0) / 0.25)        # continuidad en media celda

    def test_resistencias_negativas_o_incompletas_se_rechazan(self):
        with self.assertRaises(ValueError):
            resistencia_de_capa({"resistencia_m2_k_w": -1.0})
        with self.assertRaises(ValueError):
            resistencia_de_capa({"espesor_m": 0.1})


class RetiroTest(unittest.TestCase):
    def test_referencia_del_retiro_colocacion_o_inicio_de_obra(self):
        v = vertida(1, inicio_h=24.0)
        base = {"nombre": "moldaje", "espesor_m": 0.02, "conductividad_w_m_k": 0.1, "retiro_h": 72.0}
        desde_colocacion = proteccion_desde_config(v, {**base, "retiro_desde": "colocacion"})
        desde_obra = proteccion_desde_config(v, {**base, "retiro_desde": "inicio_obra"})
        self.assertEqual(desde_colocacion.retiro_s, (24.0 + 72.0) * H)
        self.assertEqual(desde_obra.retiro_s, 72.0 * H)
        # Valor por defecto: colocación.
        self.assertEqual(proteccion_desde_config(v, base).retiro_s, (24.0 + 72.0) * H)
        with self.assertRaises(ValueError):    # retiro antes de colocar el hormigón
            proteccion_desde_config(v, {**base, "retiro_h": 10.0, "retiro_desde": "inicio_obra"})
        with self.assertRaises(ValueError):
            proteccion_desde_config(v, {**base, "retiro_desde": "otro"})

    def test_retiro_ifc_usa_RetiroMoldaje_h_y_nunca_retira_con_null(self):
        v = vertida(1, propiedades={"RetiroMoldaje_h": 72.0})
        cfg = {"nombre": "moldaje", "resistencia_m2_k_w": 0.2}
        self.assertEqual(proteccion_desde_config(v, {**cfg, "retiro_h": "ifc"}).retiro_s, 72.0 * H)
        self.assertIsNone(proteccion_desde_config(v, {**cfg, "retiro_h": None}).retiro_s)
        with self.assertRaises(ValueError):
            proteccion_desde_config(vertida(1), {**cfg, "retiro_h": "ifc"})

    def simular(self, retiro_h, *, h_aire=1000.0, resistencia=5.0, horas=30.0):
        grilla = grilla_bloque()
        protecciones = () if retiro_h == "sin" else (Proteccion(
            1, "moldaje", CARAS, resistencia, None if retiro_h is None else retiro_h * H),)
        bordes = ModeloBordes(grilla, robin(h_aire, 10.0), protecciones=protecciones)
        resultado = simular_vertidas(
            grilla=grilla, vertidas=[vertida(1, 0.0, 40.0)], conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            hidratacion=ModeloPulsoCosenoidal(0.0, 12 * H), tiempo_final_s=horas * H,
            paso_tiempo_max_s=1800.0, tiempos_muestreo_s=[5 * H, 10 * H, 20 * H], bordes=bordes)
        return resultado

    def indice(self, resultado, horas):
        return int(np.argmin(np.abs(resultado.tiempos_s - horas * H)))

    def test_retiro_conserva_el_estado_y_luego_enfria_mas(self):
        con_retiro = self.simular(10.0)
        nunca = self.simular(None)
        n = self.indice(con_retiro, 10.0)
        # Hasta el retiro, ambas corridas son la misma (la protección sigue puesta).
        self.assertTrue(np.allclose(con_retiro.temperaturas_c[n], nunca.temperaturas_c[self.indice(nunca, 10.0)]))
        self.assertEqual(len(con_retiro.retiros), 1)
        self.assertAlmostEqual(con_retiro.retiros[0][0], 10.0 * H)
        # Sin la protección el bloque pierde más calor.
        self.assertLess(con_retiro.temperaturas_maximas_c[-1], nunca.temperaturas_maximas_c[-1])
        # El balance de energía se mantiene a través del retiro.
        self.assertLess(abs(con_retiro.error_balance_j), 1e-6 * abs(con_retiro.energia_final_j))

    def test_el_paso_estable_se_recalcula_cuando_cambian_las_conductancias(self):
        resultado = self.simular(10.0)
        _, _, antes, despues = resultado.retiros[0]
        self.assertLess(despues, antes)       # más conductancia -> paso máximo estable menor
        self.assertAlmostEqual(resultado.ultimo_paso_maximo_estable_s, despues)

    def test_bordes_sin_extras_reproducen_las_fronteras_anteriores(self):
        grilla = grilla_bloque()
        comunes = dict(grilla=grilla, vertidas=[vertida(1, 0.0, 40.0)], conductividad_w_m_k=2.0,
                       densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
                       hidratacion=ModeloPulsoCosenoidal(500.0, 12 * H), tiempo_final_s=20 * H,
                       paso_tiempo_max_s=1800.0, tiempos_muestreo_s=[5 * H, 10 * H])
        anterior = simular_vertidas(fronteras={c: robin(8.0, 20.0) for c in CARAS}, **comunes)
        nuevo = simular_vertidas(bordes=ModeloBordes(grilla, robin(8.0, 20.0)), **comunes)
        for a, b in zip(anterior.temperaturas_c, nuevo.temperaturas_c):
            self.assertTrue(np.array_equal(a, b))
        self.assertEqual(anterior.energia_frontera_j, nuevo.energia_frontera_j)

    def test_extremos_usan_la_resistencia_de_proteccion_y_el_suelo(self):
        grilla = grilla_bloque()
        suelo = ModeloSuelo(10.0, 1.5, 0.75)
        protecciones = (Proteccion(1, "moldaje", ("x+",), 1.0, None),)
        bordes = ModeloBordes(grilla, robin(10.0, 20.0), suelo=suelo, protecciones=protecciones)
        extremos = SeguimientoExtremos(grilla)
        simular_vertidas(
            grilla=grilla, vertidas=[vertida(1, 0.0, 40.0)], conductividad_w_m_k=2.0,
            densidad_kg_m3=2400.0, calor_especifico_j_kg_k=900.0,
            hidratacion=ModeloPulsoCosenoidal(0.0, 12 * H), tiempo_final_s=H, paso_tiempo_max_s=600.0,
            bordes=bordes, seguimiento=extremos)
        s = extremos.serie("global")
        # Caras a t=0 (T=40): aire sin protección, x+ protegida (R''=1) y base sobre suelo (R=0.5).
        t_aire = 40.0 + (20.0 - 40.0) * 0.25 / (0.25 + 0.1)
        t_x_mas = 40.0 + (20.0 - 40.0) * 0.25 / (0.25 + 1.0 + 0.1)
        t_base = 40.0 + (10.0 - 40.0) * 0.25 / (0.25 + 0.5)
        self.assertLess(t_aire, t_base)
        self.assertLess(t_base, t_x_mas)       # la protección deja la cara del hormigón más caliente
        self.assertAlmostEqual(s["fmin"][0], t_aire, places=6)
        self.assertAlmostEqual(s["fmax"][0], t_x_mas, places=6)
        self.assertTrue(np.isfinite(s["contraste_abs_c"][-1]))

    def test_crear_modelo_desde_config_y_errores(self):
        grilla = grilla_bloque()
        v = [vertida(1, propiedades={"RetiroMoldaje_h": 72.0})]
        cfg = {"suelo": {"activo": True, "temperatura_c": 15.0, "conductividad_w_m_k": 1.0,
                         "espesor_efectivo_m": 0.5},
               "protecciones_por_vertida": {"1": [
                   {"nombre": "moldaje", "caras": ["x+"], "resistencia_m2_k_w": 0.3, "retiro_h": "ifc"}]}}
        modelo = crear_modelo_bordes(grilla, v, cfg, robin(8.0, 20.0))
        self.assertEqual(modelo.suelo.h_equivalente_w_m2_k, 2.0)
        self.assertEqual(modelo.instantes_retiro_s(100 * H), [72.0 * H])
        with self.assertRaises(ValueError):
            crear_modelo_bordes(grilla, v, {"protecciones_por_vertida": {"9": []}}, robin(8.0, 20.0))


def serie_horaria(horas, inicio="2026-11-11T12:00:00Z", base=15.0):
    t0 = a_utc(inicio)
    tiempos = tuple(t0 + timedelta(hours=h) for h in range(horas + 1))
    valores = np.array([base + 5.0 * np.sin(2 * np.pi * h / 24.0) for h in range(horas + 1)])
    return SerieTemperatura(tiempos, valores, tuple("pronostico" for _ in tiempos))


class ClimaTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.raiz = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def snapshot(self, horas=48, momento=None):
        serie = serie_horaria(horas)
        return guardar_snapshot(
            self.raiz / "data" / "weather", bruto=b'{"crudo": true}', nombre_original="original.json",
            serie=serie, fuente="prueba", modelo="m", latitud=-30.0, longitud=-71.0,
            descargado_utc=momento or datetime(2026, 11, 11, 12, 0, tzinfo=timezone.utc)), serie

    def test_snapshot_conserva_original_serie_y_metadatos_y_no_sobrescribe(self):
        c1, serie = self.snapshot()
        c2, _ = self.snapshot()                           # misma descarga: no pisa la anterior
        self.assertNotEqual(c1, c2)
        self.assertEqual(len(listar_snapshots(self.raiz / "data" / "weather")), 2)
        meta = json.loads((c1 / "metadata.json").read_text(encoding="utf-8"))
        for clave in ("source", "model", "latitude", "longitude", "issued_at", "retrieved_at",
                      "valid_from", "valid_to", "units", "original_sha256"):
            self.assertIn(clave, meta)
        self.assertIsNone(meta["issued_at"])              # emisión no informada: se dice explícitamente
        self.assertEqual(meta["units"]["temperature_2m"], "degC")
        self.assertEqual((c1 / "original.json").read_bytes(), b'{"crudo": true}')
        leido = leer_snapshot(c1)
        self.assertTrue(np.allclose(leido.serie.temperatura_c, serie.temperatura_c, atol=1e-4))
        (c1 / "original.json").write_bytes(b"adulterado")
        with self.assertRaises(ValueError):
            leer_snapshot(c1)

    def test_cobertura_informa_el_intervalo_disponible(self):
        serie = serie_horaria(48)
        inicio = a_utc("2026-11-11T12:00:00Z")
        completa = cobertura(serie, inicio, inicio + timedelta(hours=48))
        self.assertTrue(completa.cubre_completo)
        corta = cobertura(serie, inicio, inicio + timedelta(hours=100))
        self.assertFalse(corta.cubre_completo)
        self.assertAlmostEqual(corta.horas_faltantes_final, 52.0)
        self.assertIn("INCOMPLETA", corta.resumen())
        self.assertIn("2026-11-13T12:00:00Z", corta.resumen())

    def test_extension_se_marca_y_no_altera_el_pronostico_original(self):
        serie = serie_horaria(48)
        inicio = a_utc("2026-11-11T12:00:00Z")
        extendida = extender_serie(serie, inicio + timedelta(hours=96), "ciclo_diario")
        self.assertEqual(extendida.fraccion_por_origen(), {"pronostico": 49, "extension:ciclo_diario": 48})
        self.assertTrue(np.array_equal(extendida.temperatura_c[:49], serie.temperatura_c))
        # La cola repite el último ciclo de 24 h.
        self.assertAlmostEqual(extendida.temperatura_c[49 + 5], serie.temperatura_c[25 + 5])
        constante = extender_serie(serie, inicio + timedelta(hours=60), "ultimo_valor")
        self.assertTrue(np.allclose(constante.temperatura_c[49:], serie.temperatura_c[-1]))
        with self.assertRaises(ValueError):
            extender_serie(serie_horaria(10), inicio + timedelta(hours=48), "ciclo_diario")

    def config_snapshot(self, carpeta, extension=None):
        return {"fuente": "snapshot", "inicio_obra_utc": "2026-11-11T12:00:00Z",
                "snapshot": {"carpeta": str(carpeta.relative_to(self.raiz)), "extension": extension}}

    def test_loader_con_cobertura_incompleta_exige_extension_explicita(self):
        carpeta, serie = self.snapshot(48)
        with self.assertRaises(ValueError) as error:
            construir_aire(self.config_snapshot(carpeta), tiempo_final_s=96 * H, raiz_repo=self.raiz)
        self.assertIn("INCOMPLETA", str(error.exception))
        aire = construir_aire(self.config_snapshot(carpeta, "ciclo_diario"),
                              tiempo_final_s=96 * H, raiz_repo=self.raiz)
        self.assertEqual(aire.descripcion["extension"], "ciclo_diario")
        self.assertEqual(aire.descripcion["extension_desde_utc"], "2026-11-13T12:00:00Z")
        self.assertIn("extension:ciclo_diario", aire.descripcion["puntos_por_origen"])
        self.assertAlmostEqual(aire.funcion(10 * H), serie.temperatura_c[10], places=3)
        self.assertAlmostEqual(aire.funcion(60 * H), serie.temperatura_c[12], places=3)   # cola = ciclo repetido

    def test_loader_dentro_de_la_cobertura_no_extiende_y_fuera_falla(self):
        carpeta, _ = self.snapshot(48)
        aire = construir_aire(self.config_snapshot(carpeta), tiempo_final_s=40 * H, raiz_repo=self.raiz)
        self.assertTrue(aire.cobertura.cubre_completo)
        with self.assertRaises(ValueError):
            aire.funcion(60 * H)                     # nunca inventa valores fuera de la serie

    def test_loader_sintetico_se_identifica_y_cambiar_la_serie_cambia_el_aire(self):
        cfg = {"fuente": "sintetica", "inicio_obra_utc": "2026-11-11T12:00:00Z",
               "sintetica": {"media_c": 16.0, "amplitud_diaria_c": 6.0}}
        a = construir_aire(cfg, tiempo_final_s=48 * H, raiz_repo=self.raiz)
        cfg2 = {**cfg, "sintetica": {"media_c": 26.0, "amplitud_diaria_c": 6.0}}
        b = construir_aire(cfg2, tiempo_final_s=48 * H, raiz_repo=self.raiz)
        self.assertTrue(a.es_sintetica)
        self.assertIn("SINT", a.descripcion["aviso"])
        self.assertAlmostEqual(b.funcion(7 * H) - a.funcion(7 * H), 10.0)
        constante = construir_aire(None, tiempo_final_s=H, raiz_repo=self.raiz, temperatura_constante_c=20.0)
        self.assertEqual(constante.funcion(123.0), 20.0)
        self.assertFalse(constante.es_sintetica)

    def test_normalizar_serie_valida_orden_duplicados_faltantes_y_unidades(self):
        t = ["2026-11-11T00:00", "2026-11-11T01:00", "2026-11-11T02:00"]
        self.assertEqual(len(normalizar_serie(t, [1.0, 2.0, 3.0]).tiempos_utc), 3)
        with self.assertRaises(ValueError):
            normalizar_serie(t, [1.0, None, 3.0])                       # faltante
        self.assertEqual(len(normalizar_serie(t, [1.0, None, 3.0], faltantes="descartar").tiempos_utc), 2)
        with self.assertRaises(ValueError):
            normalizar_serie([t[0], t[0], t[2]], [1.0, 2.0, 3.0])       # duplicado
        with self.assertRaises(ValueError):
            normalizar_serie([t[1], t[0], t[2]], [1.0, 2.0, 3.0])       # desordenado
        with self.assertRaises(ValueError):
            normalizar_serie(t, [1.0, 2.0, 3.0], unidad="degF")         # unidad
        kelvin = normalizar_serie(t, [273.15, 283.15, 293.15], unidad="K")
        self.assertTrue(np.allclose(kelvin.temperatura_c, [0.0, 10.0, 20.0]))

    def test_tipo_de_serie_se_registra_y_se_valida(self):
        serie = serie_horaria(30)
        for tipo in ("pronostico", "reanalisis", "observacion", "sintetica"):
            carpeta = guardar_snapshot(
                self.raiz / tipo, bruto=b"x", nombre_original="original.json", serie=serie,
                fuente="p", modelo=None, latitud=None, longitud=None, tipo=tipo,
                descargado_utc=datetime(2026, 11, 11, tzinfo=timezone.utc))
            meta = json.loads((carpeta / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["series_kind"], tipo)
            self.assertEqual(meta["synthetic"], tipo == "sintetica")
            self.assertEqual(set(leer_snapshot(carpeta).serie.origen), {tipo})
        with self.assertRaises(ValueError):
            guardar_snapshot(self.raiz / "x", bruto=b"x", nombre_original="o.json", serie=serie,
                             fuente="p", modelo=None, latitud=None, longitud=None, tipo="otro",
                             descargado_utc=datetime(2026, 11, 11, tzinfo=timezone.utc))

    def test_serie_sintetica_determinista(self):
        a = serie_sintetica("2026-11-11T12:00:00Z", 24.0)
        b = serie_sintetica("2026-11-11T12:00:00Z", 24.0)
        self.assertTrue(np.array_equal(a.temperatura_c, b.temperatura_c))
        self.assertEqual(set(a.origen), {"sintetica"})

    RESPUESTA = json.dumps({"utc_offset_seconds": 0, "hourly_units": {"time": "iso8601", "temperature_2m": "°C"},
                            "hourly": {"time": ["2026-11-11T00:00", "2026-11-11T01:00", "2026-11-11T02:00"],
                                       "temperature_2m": [14.0, 13.5, None]}}).encode("utf-8")

    def test_open_meteo_parsea_url_y_snapshot_sin_red(self):
        serie = parsear_respuesta(self.RESPUESTA)
        self.assertEqual(len(serie.tiempos_utc), 2)          # el valor nulo se descarta
        url = construir_url(-30.6, -71.2)
        self.assertIn("temperature_2m", url)
        self.assertIn("timezone=UTC", url)
        sitio = Sitio(3, "Ovalle / Ecos del Mar", "Coquimbo", -30.6, -71.2, "prueba")
        carpeta = descargar_y_guardar(sitio, self.raiz / "w", obtener=lambda _u: self.RESPUESTA,
                                      ahora_utc=datetime(2026, 11, 11, 12, tzinfo=timezone.utc))
        meta = json.loads((carpeta / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta["latitude"], -30.6)
        self.assertFalse(meta["synthetic"])
        with self.assertRaises(ValueError):                  # unidades distintas se rechazan
            parsear_respuesta(self.RESPUESTA.replace(rb"\u00b0C", b"degF"))

    def test_sin_coordenadas_se_bloquea_la_consulta_real_antes_de_la_red(self):
        sitio = sitio_desde_config({"grupo": 3, "nombre": "Ovalle / Ecos del Mar", "region": "Coquimbo"})
        self.assertFalse(sitio.tiene_coordenadas)
        llamadas = []
        with self.assertRaises(SitioSinCoordenadas) as error:
            descargar_y_guardar(sitio, self.raiz, obtener=lambda u: llamadas.append(u) or b"")
        self.assertEqual(llamadas, [])
        self.assertIn("sitios.yaml", str(error.exception))


if __name__ == "__main__":
    unittest.main()
