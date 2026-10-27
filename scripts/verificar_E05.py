"""Verificacion E05: soluciones exactas 2D/3D y fuente uniforme adiabatica.

Ejecutable con "Ejecutar archivo de Python" de VS Code (no requiere PYTHONPATH).
Edita PARAMETROS EDITABLES y vuelve a ejecutar. Termina con codigo 1 si falla
alguna verificacion (no usa assert: sigue activa con ``python -O``).

Casos (Dirichlet T0 en todas las caras, Qdot = 0 en los sinusoidales):
  2D: T = T0 + A sin(pi x/Lx) sin(pi y/Ly) exp(-alfa lambda t)
  3D: idem con sin(pi z/Lz); lambda = suma de (pi/L)^2 por eje; alfa = k/(rho cp)
  Fuente: dominio adiabatico, T = T_ini + q0 t/(rho cp), dE = q0 V_total t.
Las formulas y la tolerancia viven en mcplan.thermal.verificacion (las mismas
que usan los tests).
"""
from pathlib import Path
import sys

RAIZ_PROYECTO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ_PROYECTO / "src"))

from mcplan.thermal.verificacion import (  # noqa: E402
    verificar_fuente_uniforme, verificar_sinusoidal)

# ---------------------------------------------------------------------------
# PARAMETROS EDITABLES (unidades: m, s, °C, W/(m K), kg/m3, J/(kg K), W/m3)
# ---------------------------------------------------------------------------
CONDUCTIVIDAD_W_M_K = 2.0
DENSIDAD_KG_M3 = 2400.0
CALOR_ESPECIFICO_J_KG_K = 900.0

T0_C = 20.0                    # temperatura base y de las caras Dirichlet [°C]
AMPLITUD_C = 30.0              # A [°C]

CASO_2D = dict(forma=(30, 20), longitudes_m=(0.6, 0.4),            # celdas, [Lx, Ly] m
               tiempo_final_s=6 * 3600.0)
CASO_3D = dict(forma=(20, 16, 12), longitudes_m=(0.5, 0.4, 0.3),   # celdas, [Lx, Ly, Lz] m
               tiempo_final_s=3 * 3600.0)
FRACCION_PASO_ESTABLE = 0.9    # dt = fraccion * limite de estabilidad
CONSTANTE_TOLERANCIA = 1.0     # tol = C * (estimacion espacial + temporal con dt efectivo)

FUENTE = dict(forma=(6, 5, 4), dx_m=(0.1, 0.1, 0.1),               # celdas, m
              T_inicial_c=20.0, q0_w_m3=500.0, tiempo_final_s=6 * 3600.0)
TOLERANCIA_RELATIVA_FUENTE = 1e-10   # redondeo de maquina: la fuente es exacta en Euler


def main() -> int:
    propiedades = dict(k_w_m_k=CONDUCTIVIDAD_W_M_K, rho_kg_m3=DENSIDAD_KG_M3,
                       cp_j_kg_k=CALOR_ESPECIFICO_J_KG_K)
    resultados = []
    for caso in (CASO_2D, CASO_3D):
        resultados.append(verificar_sinusoidal(
            T0_c=T0_C, amplitud_c=AMPLITUD_C, fraccion_paso=FRACCION_PASO_ESTABLE,
            constante_tolerancia=CONSTANTE_TOLERANCIA, **propiedades, **caso))
    resultados.append(verificar_fuente_uniforme(
        tolerancia_relativa=TOLERANCIA_RELATIVA_FUENTE, **propiedades, **FUENTE))

    for r in resultados:
        d = r.detalle
        print(f"\n== {r.nombre}: {'APROBADO' if r.aprobado else 'FALLA'}")
        print(f"   error absoluto maximo T : {r.error_max_c:.4e} °C")
        print(f"   criterio (tolerancia)   : <= {r.tolerancia_c:.4e} °C")
        if "error_normalizado" in d:
            en = d["error_normalizado"]
            print(f"   error / |A|             : {'n/a (A = 0)' if en is None else format(en, '.3e')}")
            print(f"   estimacion espacial/temporal : {d['estimacion_espacial_c']:.3e} / "
                  f"{d['estimacion_temporal_c']:.3e} °C (criterio = suma de magnitudes)")
            print(f"   error / tolerancia      : {d['error_sobre_tolerancia']:.3f}")
            print(f"   decaimiento exp(-a*l*t) : {d['factor_decaimiento']:.4f}  "
                  f"(T_max exacta {d['T_max_exacta_c']:.3f} °C, numerica {d['T_max_numerica_c']:.3f} °C)")
            print(f"   dt = {d['dt_s']:.2f} s (limite {d['dt_limite_s']:.2f} s), {d['pasos']} pasos")
        else:
            print(f"   T esperada              : {d['T_esperada_c']:.6f} °C")
            e_ok = d["error_energia_relativo"] <= d["tol_energia_relativa"]
            print(f"   dE numerica / esperada  : {d['dE_j']:.6e} / {d['dE_esperado_j']:.6e} J "
                  f"(q0*V*t, V = {d['V_total_m3']:.3f} m3)")
            print(f"   error absoluto energia  : {d['error_energia_j']:.2e} J "
                  f"(escala {d['denominador_energia_j']:.3e} J)")
            print(f"   error energia/escala    : {d['error_energia_relativo']:.2e} "
                  f"(criterio <= {d['tol_energia_relativa']:.0e}) -> {'OK' if e_ok else 'FALLA'}")
    aprobados = all(r.aprobado for r in resultados)
    aprobados = aprobados and resultados[-1].detalle["error_energia_relativo"] <= \
        resultados[-1].detalle["tol_energia_relativa"]
    print("\nRESULTADO GLOBAL:", "TODAS LAS VERIFICACIONES APROBADAS" if aprobados else "HAY VERIFICACIONES FALLIDAS")
    return 0 if aprobados else 1


if __name__ == "__main__":
    raise SystemExit(main())
