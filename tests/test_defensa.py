"""Pruebas de la búsqueda de la defensa: horarios, D, horizonte limitado y métricas auxiliares."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from mcplan.geometry.voxelizacion import GrillaVoxel
from mcplan.thermal.decision import SeguimientoMetricasDecision
from mcplan.thermal.diffusion_nd import DifusionCartesiana, dirichlet

import defensa_busqueda as db


class PlanDefensaTest(unittest.TestCase):
    def test_vertida_continua_no_tiene_esperas_y_D_es_volumen_sobre_tasa(self):
        h = db.construir_plan({"id": "u", "tasa_m3_h": 36.0})
        self.assertAlmostEqual(h["D_h"], sum(db.VOLUMENES_M3.values()) / 36.0)
        for i in range(1, 5):
            self.assertAlmostEqual(h["activaciones_h"][i] + h["duraciones_h"][i], h["activaciones_h"][i + 1])

    def test_D_incluye_esperas_y_postergacion_inicial(self):
        h = db.construir_plan({"id": "w", "tasa_m3_h": 36.0, "espera_h": [2.0, 0.0, 3.0, 1.0], "inicio_primera_h": 12.0})
        self.assertAlmostEqual(h["D_h"], 12.0 + sum(db.VOLUMENES_M3.values()) / 36.0 + 6.0)

    def test_rechaza_solape_con_un_solo_frente(self):
        with self.assertRaises(ValueError):
            db.construir_plan({"id": "x", "activaciones_h": {1: 0, 2: 1, 3: 9, 4: 12, 5: 15},
                               "duraciones_h": {1: 4, 2: 4, 3: 2, 4: 2, 5: 2}})


class ObservadorHorizonteTest(unittest.TestCase):
    def setUp(self):
        dx = 0.1
        forma = (4, 1, 1)
        self.grilla = GrillaVoxel((0.0, 0.0, 0.0), (dx, dx, dx), np.ones(forma, dtype=np.int32), {1: 4 * dx**3})
        self.solver = DifusionCartesiana(
            forma=forma, dx_m=(dx, dx, dx), conductividad_w_m_k=2.0, densidad_kg_m3=2400.0,
            calor_especifico_j_kg_k=900.0,
            fronteras={c: dirichlet(0.0) for c in ("x-", "x+", "y-", "y+", "z-", "z+")},
            mascara_activa=np.ones(forma, dtype=bool))

    def test_tiempo_max_ignora_estados_posteriores(self):
        corto = SeguimientoMetricasDecision(self.grilla, 0.025, tiempo_max_s=3600.0)
        largo = SeguimientoMetricasDecision(self.grilla, 0.025)
        for t, escala in ((0.0, 1.0), (3600.0, 1.0), (7200.0, 3.0)):
            campo = (escala * np.array([20.0, 30.0, 40.0, 50.0])).reshape(4, 1, 1)
            corto.registrar(t, campo, self.solver)
            largo.registrar(t, campo, self.solver)
        self.assertAlmostEqual(corto.resumen()["M_delta_c"], 40.0)
        self.assertGreater(largo.resumen()["M_delta_c"], 100.0)
        self.assertEqual(corto.resumen()["n_muestras"], 2)

    def test_metricas_auxiliares_tienen_nombre_propio_y_no_sustituyen_a_S25(self):
        obs = SeguimientoMetricasDecision(self.grilla, 0.025)
        obs.registrar(0.0, np.array([20.0, 30.0, 40.0, 50.0]).reshape(4, 1, 1), self.solver)
        r = obs.resumen()
        self.assertAlmostEqual(r["M_delta_c"], 40.0)             # S25: 50 - 10
        self.assertAlmostEqual(r["M_delta_cara_c"], 50.0)        # cara a 0 °C
        self.assertAlmostEqual(r["M_delta_centro_voxel_c"], 30.0)  # 50 - 20 (centro del voxel de borde)


class PropuestaBTest(unittest.TestCase):
    """La configuración documentada de la propuesta B coincide con lo que ejecuta el código."""

    def test_config_B_coincide_con_el_calendario_y_el_plan_ejecutados(self):
        import json
        import defensa_comprobaciones as dc
        raiz = Path(__file__).resolve().parents[1]
        cfg = json.loads((raiz / "configs/defensa/propuesta_B.json").read_text(encoding="utf-8"))
        self.assertEqual(cfg["protecciones"]["manta_PUR"]["retiro_capas_dias_desde_t0"], [56, 77, 98, 119])
        tareas = [t for t in dc.tareas() if t["grupo"] == "retB"]
        self.assertEqual(len(tareas), 48)
        spec = tareas[0]["spec"]
        self.assertEqual([h / 24 for h in spec["manta_capas"]["retiros_h"]], [56, 77, 98, 119])
        self.assertAlmostEqual(spec["manta_capas"]["espesor_capa_m"] * spec["manta_capas"]["n"],
                               cfg["protecciones"]["manta_PUR"]["espesor_inicial_m"])
        self.assertAlmostEqual(spec["aislacion_inferior"]["espesor_m"], cfg["protecciones"]["aislacion_inferior"]["espesor_m"])
        self.assertEqual(spec["temperatura_colocacion_c"], cfg["temperatura_colocacion_c"])
        horarios = db.construir_plan(spec)
        self.assertAlmostEqual(spec["tasa_m3_h"], cfg["colocacion"]["tasa_m3_h"])
        self.assertAlmostEqual(horarios["D_h"], cfg["colocacion"]["D_h"], places=2)

    def test_tiempos_frente_respetan_inicio_fin_y_caudal_constante(self):
        import defensa_comprobaciones as dc
        h = db.construir_plan({"id": "t", "tasa_m3_h": 36.0})

        class V:  # vertida mínima
            def __init__(self, i):
                self.vertida_id = i
        ids = np.zeros((4, 4, 4), dtype=np.int32)
        ids[:, :, :2] = 1
        ids[:, :, 2:] = 2
        grilla = type("G", (), {"forma": ids.shape, "vertida_id": ids})()
        t = dc.tiempos_frente(grilla, [V(1), V(2)], h)
        n1 = int((ids == 1).sum())
        self.assertGreaterEqual(t[ids == 1].min(), h["activaciones_h"][1] * 3600.0 - 1e-6)
        self.assertLessEqual(t[ids == 1].min(), h["activaciones_h"][1] * 3600.0 + h["duraciones_h"][1] * 3600.0 / n1 + 300.0)
        self.assertLessEqual(t[ids == 1].max(), (h["activaciones_h"][1] + h["duraciones_h"][1]) * 3600.0 + 1e-6)
        self.assertTrue(np.all(np.diff(np.sort(t[ids == 1])) >= 0))
        self.assertGreaterEqual(t[ids == 2].min(), h["activaciones_h"][2] * 3600.0 - 1e-6)


if __name__ == "__main__":
    unittest.main()
