"""Genera los archivos de planes (JSON) del cribado y de la comparación final.

    python scripts/defensa_generar_planes.py cribado
    python scripts/defensa_generar_planes.py tiempos <id_paquete>
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys

RAIZ = Path(__file__).resolve().parents[1]
DESTINO = RAIZ / "configs/defensa"
VOL_TOTAL = 636.776
R_ALTA = 63.7   # caudal medio implícito en la secuencia de 10 h del compañero
R_EQUIPO = 36.0  # supuesto de configs/evaluacion_final.jsonc (una cuadrilla)

SUELOS = {"L2m": {"T": 20.0, "k": 1.5, "L": 2.0}, "L0_8m": {"T": 20.0, "k": 1.5, "L": 0.8}}


def paquete(pid, *, tcol, manta_mm, aisl_mm=0, suelo="L2m", retiro_manta=None, manta_inicio=15.0):
    p = {"tcol": tcol, "manta_mm": manta_mm, "aisl_mm": aisl_mm, "suelo": suelo,
         "retiro_manta": retiro_manta, "manta_inicio": manta_inicio}
    p["id"] = pid
    return p


def spec_plan(pid, cat, pq, *, tasa=R_ALTA, espera=0.0, inicio=0.0, **extra):
    s = {"id": pid, "categoria": cat, "tasa_m3_h": tasa, "espera_h": espera, "inicio_primera_h": inicio,
         "temperatura_colocacion_c": pq["tcol"],
         "moldaje": {"vertidas": [1, 2, 5], "espesor_m": 0.018, "k": 0.1},
         "suelo": SUELOS[pq["suelo"]]}
    if pq["manta_mm"]:
        s["manta"] = {"espesor_m": pq["manta_mm"] / 1000.0, "k": 0.025, "inicio_h": pq["manta_inicio"],
                      "retiro_h": pq["retiro_manta"]}
    if pq["aisl_mm"]:
        s["aislacion_inferior"] = {"espesor_m": pq["aisl_mm"] / 1000.0, "k": 0.04}
    s.update(extra)
    return s


def cribado():
    planes = []
    for suelo, tcol, mm, aisl in itertools.product(SUELOS, (10.0, 15.0, 18.0), (40, 60, 80), (0, 50)):
        pid = f"c_{suelo}_T{tcol:.0f}_PUR{mm}_AI{aisl}"
        planes.append(spec_plan(pid, "cribado", paquete(pid, tcol=tcol, manta_mm=mm, aisl_mm=aisl, suelo=suelo)))
    return planes


def cribado2():
    """Segunda ronda: más aislación, para localizar el máximo tardío de ΔT (horizonte 56 d)."""
    planes = []
    for suelo, tcol, mm, aisl in itertools.product(SUELOS, (10.0, 15.0), (60, 80, 100), (0, 50, 100)):
        if aisl == 0 and (suelo != "L2m" or mm == 60):
            continue
        pid = f"d_{suelo}_T{tcol:.0f}_PUR{mm}_AI{aisl}"
        planes.append(spec_plan(pid, "cribado2", paquete(pid, tcol=tcol, manta_mm=mm, aisl_mm=aisl, suelo=suelo)))
    return planes


def secuencia_fija(paso_h, dur_h):
    """Cinco vertidas, una cada ``paso_h`` con duración ``dur_h`` (una cuadrilla, sin traslape)."""
    return {"activaciones_h": {str(i): (i - 1) * paso_h for i in range(1, 6)},
            "duraciones_h": {str(i): dur_h for i in range(1, 6)}}


def final():
    """Comparación final: paquetes P2/P3, secuencias y referencias del equipo."""
    P2 = paquete("P2", tcol=10.0, manta_mm=80, aisl_mm=50)
    P3 = paquete("P3", tcol=15.0, manta_mm=100, aisl_mm=100)
    P1 = paquete("P1", tcol=10.0, manta_mm=40, aisl_mm=0)
    planes = [
        spec_plan("f_P2_L2m_U63", "vertida_unica", P2, tasa=R_ALTA,
                  descripcion="Colocación continua de las cinco partes a 63,7 m³/h, sin esperas."),
        spec_plan("f_P2_L2m_U36", "mas_rapido_equipo", P2, tasa=R_EQUIPO,
                  descripcion="Colocación continua a 36 m³/h (supuesto de una cuadrilla del equipo)."),
        spec_plan("f_P2_L2m_4_5h", "mas_rapido", P2, **secuencia_fija(4.5, 4.5),
                  descripcion="Una vertida cada 4,5 h (más rápida que la propuesta de 6 h)."),
        spec_plan("f_P2_L2m_6h", "propuesto", P2, **secuencia_fija(6.0, 6.0),
                  descripcion="Plan propuesto del equipo (IFC): una vertida cada 6 h de 6 h."),
        spec_plan("f_P2_L2m_8h", "mas_lento", P2, **secuencia_fija(8.0, 8.0),
                  descripcion="Una vertida cada 8 h."),
        spec_plan("f_P2_L2m_12h", "mas_lento", P2, **secuencia_fija(12.0, 6.0),
                  descripcion="Una vertida cada 12 h (6 h de colocación + 6 h de espera)."),
        spec_plan("f_P2_L0_8m_U63", "vertida_unica", dict(P2, suelo="L0_8m"), tasa=R_ALTA,
                  descripcion="Igual a f_P2_L2m_U63 con suelo efectivo de 0,8 m (E10)."),
        spec_plan("f_P2_L0_8m_6h", "propuesto", dict(P2, suelo="L0_8m"), **secuencia_fija(6.0, 6.0),
                  descripcion="Plan de 6 h con suelo efectivo de 0,8 m."),
        spec_plan("f_P3_L2m_U63", "vertida_unica", P3, tasa=R_ALTA,
                  descripcion="Colocación continua con T de colocación 15 °C, PUR 100 mm y aislación inferior 100 mm."),
        spec_plan("f_P3_L0_8m_U63", "vertida_unica", dict(P3, suelo="L0_8m"), tasa=R_ALTA,
                  descripcion="Igual con suelo efectivo de 0,8 m."),
        spec_plan("f_P1_L2m_U63", "referencia_equipo", P1, tasa=R_ALTA,
                  descripcion="Referencia del compañero: T 10 °C, PUR 40 mm, sin aislación inferior (sin retiro de manta)."),
        {"id": "f_ref_protegido_6h", "categoria": "referencia_equipo", "inicio_primera_h": 0.0,
         **secuencia_fija(6.0, 6.0), "temperatura_colocacion_c": 18.0,
         "moldaje": {"vertidas": [1, 2, 3, 4, 5], "espesor_m": 0.018, "k": 0.13},
         "manta": {"espesor_m": 0.075, "k": 0.04, "inicio_h": 0.0, "retiro_h": 360.0},
         "aislacion_inferior": {"espesor_m": 0.10, "k": 0.04},
         "suelo": {"T": 18.0, "k": 1.5, "L": 0.8},
         "descripcion": "Candidato 'propuesto_6h_protegido' del compañero (configs/verificacion_candidato_proteccion.jsonc), ahora con 8 curvas y 28/56 d."},
    ]
    return planes


def sens():
    """Sensibilidades de supuestos físicos y constructivos sobre el candidato f_P2_L2m_U63."""
    import copy
    base = spec_plan("f_P2_L2m_U63", "vertida_unica", paquete("P2", tcol=10.0, manta_mm=80, aisl_mm=50), tasa=R_ALTA)
    out = []

    def var(sufijo, fn, cat="sensibilidad"):
        s = copy.deepcopy(base)
        s["id"] = f"s_{sufijo}"
        s["categoria"] = cat
        fn(s)
        out.append(s)

    for h in (5.0, 12.0):
        var(f"h_aire_{h:g}", lambda s, h=h: s.update(h_aire=h))
    for L in (0.8, 4.0):
        var(f"suelo_L{str(L).replace('.', '_')}", lambda s, L=L: s["suelo"].update(L=L))
    for T in (16.0, 24.0):
        var(f"suelo_T{T:.0f}", lambda s, T=T: s["suelo"].update(T=T))
    for k in (0.022, 0.035):
        var(f"k_manta_{str(k).replace('.', '_')}", lambda s, k=k: s["manta"].update(k=k))
    for t in (8.0, 12.0):
        var(f"tcol_{t:.0f}", lambda s, t=t: s.update(temperatura_colocacion_c=t))
    for i in (6.0, 12.0):
        var(f"inicio_{i:g}h", lambda s, i=i: s.update(inicio_primera_h=i))
    for mi in (6.0, 24.0):
        var(f"manta_inicio_{mi:g}h", lambda s, mi=mi: s["manta"].update(inicio_h=mi))
    for r in (336.0, 672.0):
        var(f"manta_retiro_{r:g}h", lambda s, r=r: s["manta"].update(retiro_h=r), cat="control_retiro")
    var("moldaje_retiro_168h", lambda s: s["moldaje"].update(retiro_h=168.0))
    var("sin_moldaje", lambda s: s.update(moldaje=None))
    # Sensibilidad numérica del candidato y de dos competidores (más ligeros).
    comp = {
        "f_P2_L2m_U63": base,
        "d_L2m_T10_PUR60_AI50": spec_plan("d_L2m_T10_PUR60_AI50", "competidor", paquete("x", tcol=10.0, manta_mm=60, aisl_mm=50), tasa=R_ALTA),
        "f_P2_L0_8m_U63": spec_plan("f_P2_L0_8m_U63", "candidato_suelo_0_8", paquete("x", tcol=10.0, manta_mm=80, aisl_mm=50, suelo="L0_8m"), tasa=R_ALTA),
        "d_L2m_T10_PUR80_AI0": spec_plan("d_L2m_T10_PUR80_AI0", "competidor", paquete("x", tcol=10.0, manta_mm=80, aisl_mm=0), tasa=R_ALTA),
    }
    num = []
    for pid, sp in comp.items():
        for tag, kw in (("dx0_5", {"dx_m": 0.5}), ("dx0_125", {"dx_m": 0.125}),
                        ("dt300", {"dt_max_s": 300.0}), ("dt150", {"dt_max_s": 150.0})):
            s = copy.deepcopy(sp)
            s["id"] = f"n_{pid}_{tag}"
            s["categoria"] = "numerica"
            s.update(kw)
            num.append(s)
    return out, num


def main() -> int:
    DESTINO.mkdir(parents=True, exist_ok=True)
    modo = sys.argv[1]
    if modo == "cribado":
        salida = {"descripcion": "Cribado de controles (T de colocación, manta PUR sin retiro, aislación inferior, suelo) con vertida continua a 63,7 m³/h.",
                  "planes": cribado()}
        (DESTINO / "planes_cribado.json").write_text(json.dumps(salida, indent=1), encoding="utf-8")
        print(f"{len(salida['planes'])} planes -> configs/defensa/planes_cribado.json")
    elif modo == "sens":
        v, n = sens()
        (DESTINO / "planes_sens.json").write_text(json.dumps({"planes": v}, indent=1), encoding="utf-8")
        (DESTINO / "planes_numerica.json").write_text(json.dumps({"planes": n}, indent=1), encoding="utf-8")
        print(f"{len(v)} planes de sensibilidad física, {len(n)} numéricos")
    elif modo == "final":
        salida = {"descripcion": "Comparación final de planes (paquetes P1-P3 y referencias).", "planes": final()}
        (DESTINO / "planes_final.json").write_text(json.dumps(salida, indent=1), encoding="utf-8")
        print(f"{len(salida['planes'])} planes -> configs/defensa/planes_final.json")
    elif modo == "cribado2":
        salida = {"descripcion": "Segunda ronda de cribado, horizonte 56 d (continuación climática etiquetada tras 28 d).",
                  "planes": cribado2()}
        (DESTINO / "planes_cribado2.json").write_text(json.dumps(salida, indent=1), encoding="utf-8")
        print(f"{len(salida['planes'])} planes -> configs/defensa/planes_cribado2.json")
    else:
        raise SystemExit("modo desconocido")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
