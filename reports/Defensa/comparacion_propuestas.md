# Comparación de las dos propuestas (A: compañero · B: alternativa adicional)

**A** = propuesta de Vicente (commit `c1f5e26`, `reports/Semana Final/`, `configs/defensa_final_planes_100m3h.json`). **B** = alternativa adicional de este trabajo (`reports/Defensa/`, `configs/defensa/propuesta_B.json`). Ambas usan el mismo solver, IFC, curvas Qdot, pesos y trayectorias GEFS del 16-11-2026, t0 = jueves 19-11-2026 10:00 (Ovalle, UTC−3), límites 70 °C / 20 °C / P_f ≤ 5 % y S25 a 25 mm. Se conserva A tal como fue entregado.

## 1. Configuraciones (lo que cambia entre ambas)
| | **A (compañero)** | **B (alternativa)** |
|---|---|---|
| Colocación | programa de 10 h (V1 0–2, V2 2–4, V3 4–7, V4 7–9, V5 9–10 h) o continua a 100 m³/h (6,37 h) | continua progresiva a **36 m³/h**, una cuadrilla (D = 17,69 h) |
| Capacidad de colocación supuesta | hasta 100 m³/h | 36 m³/h (supuesto de una cuadrilla del equipo; sin validar) |
| Activación térmica | instantánea al inicio de cada vertida | progresiva por capas a caudal constante (comparada con instantánea: ΔT ±0,3 °C) |
| Temperatura de colocación | 10 °C | 10 °C |
| Manta PUR | 60 mm (instalada tras la colocación), sin retiro dentro del horizonte | 80 mm en 4 capas de 20 mm desde las 15 h; **retiro de 20 mm los días 56, 77, 98 y 119** |
| Aislación inferior | no | **100 mm** (propuesta térmica condicionada) |
| Suelo efectivo | 2 m, 20 °C, k = 1,5 | 2 m **y** 0,8 m (mismo 20 °C, k = 1,5) |
| Horizonte evaluado | 14 d (336 h) | 14 d, 28 d, 56 d y 141,7 d (3400 h) |
| Escenarios | 24 (8 curvas × 3 climas) | 24 (8 curvas × 3 climas), 3 interpretaciones de pesos |
| Malla / paso | 0,25 m / 600 s | 0,25 m / 600 s; escenario crítico verificado en 0,25×0,25×0,05 y 0,125×0,125×0,05 |

## 2. Resultados con métricas equivalentes (24 escenarios, interpretación de 8 curvas)
| Métrica | A: 10 h, PUR 60 mm (14 d) | B, suelo 2 m | B, suelo 0,8 m |
|---|---:|---:|---:|
| D [h] | 10,0 (6,37 a 100 m³/h) | 17,69 | 17,69 |
| Tmax peor caso [°C] | 56,30 | 56,61 | 56,58 |
| ΔT peor caso a 14 d [°C] | 18,96 | 16,41 | 16,41 |
| P_f a 14 d | 0 | 0 | 0 |
| ΔT peor caso a 56 d [°C] | no evaluado (24 esc.) | 18,01 | 18,15 |
| ΔT peor caso a 141,7 d [°C] | no evaluado (24 esc.) | 18,29 | 18,17 |
| P_f a 56 d / 141,7 d | no evaluado (24 esc.) | 0 / 0 | 0 / 0 |
| Instante del ΔT peor | 336 h (fin del horizonte) | día 60,8 (primer retiro) | día 59,8 (primer retiro) |
Para **A** el máximo de ΔT cae en la última hora simulada (336 h); más allá del día 14 no hay matriz de 24 escenarios. Como referencia **de un solo escenario** con este motor (M4_SS fría, suelo 2 m, plan análogo continuo a 63,7 m³/h, 10 °C, PUR 60 mm, sin aislación inferior ni retiro; `cribado1_resultados_por_escenario.csv`), el ΔT a 14 d es 18,96 °C (coincide con el 18,964 °C de A) y **a 28 d es 21,22 °C (> 20 °C)**: indica que A no queda demostrada más allá del día 14; no es una probabilidad conjunta. Con PUR 40 mm (primer informe de A: P_f = 0,169 a 14 d, que este motor reproduce), la matriz de 24 escenarios de un plan análogo (`f_P1_L2m_U63`, `comparacion_planes.csv`) da P_f = 0,415 a 56 d.
**Qué cambia entre A y B que impide comparar sin reservas:** recursos (100 vs 36 m³/h), horizonte (14 d vs 141,7 d), suelos (2 m vs 2 m y 0,8 m), y presencia de aislación inferior y de retiro por capas. A es más rápida (D menor) pero está demostrada solo a 14 d; B es más lenta, requiere aislación inferior y un calendario de retiro, y se evaluó hasta que el pico quedó dentro del horizonte.

## 3. Regla de retiro anterior que falla (se conserva como evidencia)
Regla **A-14d** (capas de 20 mm retiradas cada 14 d desde el día 56: días 56, 70, 84, 98 — 14-01, 28-01, 11-02, 25-02-2027), con el mismo candidato que B (AI 100 mm): ΔT máx **20,71 °C** con suelo 2 m (día 99,8; P_f = 0,081 con 4 SS y 0,040 con 8 curvas, un escenario falla) y 19,56 °C con suelo 0,8 m (margen 0,44 °C). B espacia los retiros 21 d y baja el máximo a 18,29 / 18,17 °C. Fuentes: Anexos 7 y 8 del informe, `ai100_retiro_*.csv`, `ret2_reglas_M4SS_calida.csv`, figura `propuestaB_series_A_vs_B.png`. (No confundir con la propuesta A del compañero.)

## 4. Alcance y limitaciones
- **Cumplimiento = escenarios y supuestos evaluados.** P_f = 0 se refiere al conjunto discreto de 24 escenarios (3 trayectorias GEFS × 8 curvas) con los pesos declarados; no es una probabilidad real de fisuración ni una garantía de ausencia de fisuras.
- **Climáticas:** el pronóstico GEFS solo es útil hasta ~14 d; los días 14–28 son trayectorias sin habilidad y **desde la hora 672 se repite el último ciclo diario (supuesto)**. Los retiros (días 56–119) y los resultados de 56–141,7 d dependen de ese supuesto; solo lo evaluado a 14 d está respaldado por pronóstico.
- **Numéricas:** malla 0,25 m (escalera, volúmenes por vertida con errores de ±10 % que cambian con la malla); dt 600 s sobreestima ΔT en 0,008 °C; en el escenario crítico el refinamiento de malla subió el ΔT +0,3 a +0,5 °C (observado, no aplicado a los demás escenarios ni una cota); la verificación fina cubre el pico del primer retiro con suelo 2 m y el pico base con suelo 0,8 m, no los demás escenarios.
- **Constructivas pendientes (no validadas):** material y viabilidad estructural de la aislación de 100 mm bajo la fundación; enfriamiento del hormigón a 10 °C; caudal real de 36 m³/h; producto y k de la manta y su retiro por capas; **el calendario de retiro es solo evaluado: la retirada real requerirá seguimiento térmico y un criterio de autorización basado en mediciones que no se ha validado**; fundamento del suelo efectivo (2 m / 0,8 m); significado de las 8 curvas (réplicas o mezclas).
- B **no** es una solución de obra validada.
