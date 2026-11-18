# Batería de 11 casos (264 combinaciones): plan del compañero, controles, sensibilidad a la temperatura de colocación y caso 11

Modelo `mcplan.thermal` (voxeles 0,25 m, dt máx 600 s), suelo efectivo de 2 m a 20 °C, 8 curvas de hidratación (1/8 cada una) × 3 trayectorias GEFS (pesos 0,354838709677 / 0,322580645161 / 0,322580645161, tomados del archivo; suman 1 − 1·10⁻¹²), 24 escenarios por caso, 336 h desde el 19-11-2026 10:00 local. Se reutilizan los resultados publicados cuando la configuración coincide exactamente; se recalcularon 3 escenarios por caso reutilizado y reprodujeron las métricas con diferencia 0,0 °C. Scripts: `scripts/bateria7_casos.py`, `bateria7_analisis.py`, `bateria7_temperatura.py`, `bateria7_figuras.py`; configuración `configs/defensa/bateria7_casos.json`.

**Convención de colocación.** Cada vertida se activa completa, de golpe, en su hora de inicio, con su temperatura de colocación; la edad de hidratación parte de ese instante. Los casos 2 y 5 activan las cinco vertidas en t=0: es una idealización no ejecutable (D constructiva no definida; cota por capacidad 636,776 m³ / 100 m³/h = 6,368 h). Los demás casos son secuencias.

**Manta y moldaje.** Manta PUR 60 mm instalada D + 5 h tras el fin de la colocación (15 h, 11,368 h, 29,0 h) y retirada 360 h después de instalarla (fuera de las 336 h). Moldaje 18 mm solo en V1, V2 y V5 (x±, y±), retirado 72 h después de colocada cada vertida.

## 1. Siete casos
Tabla en `resumen_bateria7.md`; datos en `tablas/resultados_por_escenario.csv`.

| Caso | P_f | M_T máx [°C] | M_Δ máx [°C] (h) |
|---|---|---|---|
| 1 Base original | 100 % | 63,27 | 45,89 (164) |
| 2 Base simultáneo | 100 % | 63,27 | 46,12 (164) |
| 3 Base largo 24 h | 100 % | 63,25 | 45,85 (164) |
| 4 Plan 1 original (10 h, 10 °C, PUR 60 mm) | 0 % | 56,30 | 18,96 (336) |
| 5 Plan simultáneo (idealizado) | 0 % | 56,12 | 18,83 (336) |
| 6 Plan largo 24 h | 0 % | 56,57 | 19,27 (336) |
| 7 Control sin manta | 100 % | 55,74 | 38,66 (164) |

La comparación 4 contra 7 aísla el efecto de la manta (mismo programa, temperatura y moldaje): M_Δ baja de 38,66 a 18,96 °C y Tmax casi no cambia. El programa de colocación (casos 1–3 y 4–6) mueve M_Δ menos de 0,5 °C.

## 2. Variantes del plan 1 con temperatura de colocación uniforme de 11, 12 y 13 °C
Solo cambia la temperatura de colocación; todo lo demás es idéntico al caso 4 (10 °C). Datos: `tablas/sensibilidad_temperatura_plan1.csv`, `tablas/plan1_T10_13_por_escenario.csv`. Figuras: `figuras/sensibilidad_T_colocacion.png`, `figuras/curvas_T10_13_*.png`.

| T colocación | M_T máx [°C] | M_Δ máx [°C] | Margen M_Δ [°C] | Escenarios con M_Δ > 20 | P_T | P_Δ | P_f conjunta |
|---|---|---|---|---|---|---|---|
| 10 | 56,298 | 18,964 | +1,036 | 0 de 24 | 0 | 0 | 0 |
| 11 | 57,255 | 19,511 | +0,489 | 0 de 24 | 0 | 0 | 0 |
| 12 | 58,211 | 20,059 | −0,059 | 1 (M4_SS frío) | 0 | 4,44 % | 4,44 % |
| 13 | 59,168 | 20,607 | −0,607 | 3 (M4_SS frío, central, cálido) | 0 | 12,50 % | 12,50 % |

- Escenario crítico: M4_SS con clima frío para M_Δ, a las 336 h, en las cuatro temperaturas; M4_SS con clima central para M_T (212–214 h). No cambia con la temperatura.
- Cambio por grado adicional (10 → 13 °C): M_T +0,96 °C/°C, M_Δ +0,55 °C/°C, prácticamente constante; por escenario, M_Δ sube exactamente 0,548 °C/°C.
- Primer valor ensayado que incumple: **12 °C**, por 0,059 °C en un escenario (M_Δ = 20,059 °C). Con P_f conjunta = 4,44 % el criterio agregado P_f ≤ 5 % todavía se cumple a 12 °C, pero el máximo de un escenario ya supera 20 °C. A 13 °C, P_f = 12,5 % > 5 %. Esto no define una temperatura límite exacta: solo se ensayaron 10, 11, 12 y 13 °C, con la malla de 0,25 m (en el escenario crítico el refinamiento de malla subió M_Δ entre 0,3 y 0,5 °C, no aplicado aquí).
- **Cumplimiento evaluado solo hasta 14 días; pico no demostrado.** En las cuatro variantes el máximo de M_Δ cae a las 336 h y la curva sigue creciendo (≈ +0,025 °C/h en el escenario crítico en las últimas 24 h). No se amplió el horizonte.
- **Alcance de la sensibilidad.** Esta sensibilidad corresponde al modelo de hidratación usado: la fuente q̇ depende solo de la edad del hormigón (curvas tabuladas con la edad medida desde la colocación), de modo que no incorpora aceleración de la hidratación por mayor temperatura (maduración/Arrhenius). El efecto de la temperatura de colocación entra solo como temperatura inicial. En un hormigón real, un colocado más caliente liberaría calor más rápido y el efecto sería probablemente mayor; esta tabla no lo cuantifica.

## 3. Caso 11: base original IFC a 10 °C
Configuración idéntica al caso 1 (inicios 0, 6, 12, 18 y 24 h, moldaje, sin manta, suelo, propiedades, hidratación, clima, malla, paso temporal); solo cambia la temperatura de colocación de las cinco vertidas, de 18 a 10 °C. 24 corridas nuevas, 336 h. Datos: `tablas/caso11_comparacion.csv`, `tablas/caso11_por_escenario.csv`; figuras: `figuras/curvas_caso11_*.png`, `figuras/3d_caso11_M4_SS_frio_164h.png` y `_336h.png`.

| Caso | Tmax máx [°C] (h) | Margen Tmax | ΔT máx [°C] (h) | Margen ΔT | Escenarios con ΔT > 20 | P_f |
|---|---|---|---|---|---|---|
| 1 Base original (18 °C) | 63,27 (192,5) | +6,73 | 45,89 (164) | −25,89 | 24 de 24 | 100 % |
| 11 Base original a 10 °C | 55,92 (202,5) | +14,08 | 38,65 (164) | −18,65 | 24 de 24 | 100 % |
| 7 Control sin manta (10 h, 10 °C) | 55,74 (194,3) | +14,26 | 38,66 (164) | −18,66 | 24 de 24 | 100 % |
| 4 Plan 1 (10 h, 10 °C, PUR 60 mm) | 56,30 (214,0) | +13,70 | 18,96 (336) | +1,04 | 0 | 0 % |

En los cuatro casos el escenario crítico es M4_SS con clima frío para ΔT y M4_SS con clima central para Tmax. Ubicaciones: Tmax cerca de (−0,13; −0,13; 1,88 m) en los casos 1, 11 y 7 y (−0,13; −0,13; 2,38 m) en el 4; ΔT crítico en la superficie S25 en (−8,98; −3,63; 1,12 m) en los casos 1, 11 y 7 y (−8,88; −3,63; 0,02 m) en el 4. El ΔT mínimo entre los 24 escenarios es 27,59 (caso 1), 20,39 (caso 11), 20,19 (caso 7) y 9,39 °C (caso 4): ningún escenario sin manta cumple.

Efectos (diferencias de máximos y, entre paréntesis, rango por escenario):
- **Temperatura de colocación, caso 1 → caso 11 (18 → 10 °C):** Tmax máx −7,35 °C (−7,48 a −7,35 por escenario) y ΔT máx −7,24 °C (−7,24 a −7,20). Es el efecto exclusivo de la temperatura (8 °C menos de colocación bajan ΔT unos 0,9 °C por grado, por encima del ~0,55 °C/°C del plan con manta).
- **Secuencia, caso 11 → caso 7 (ambos a 10 °C, sin manta):** ΔT máx +0,02 °C (−0,22 a +0,02 por escenario) y Tmax máx −0,18 °C (−0,65 a −0,15). La secuencia IFC de 24 h y el programa de 10 h son prácticamente equivalentes en este modelo.
- **Manta, caso 7 → caso 4:** ΔT máx −19,70 °C (−19,70 a −10,80 por escenario) y Tmax máx +0,56 °C. Es el efecto exclusivo de la manta, y es lo que pasa el criterio de 20 °C.

**Horizonte:** en los casos 1, 11 y 7 el máximo de ΔT cae a las 164 h (0 de 24 escenarios con el pico en las últimas 24 h): el pico queda dentro de las 336 h y la curva ya decrece. En el caso 4 el máximo cae a las 336 h en los 24 escenarios y la curva sigue creciendo: **cumplimiento evaluado hasta 14 días; pico no demostrado**. No se ejecutaron extensiones.

## 4. Salidas
`tablas/`: `resultados_por_escenario.csv` (casos 1–7), `caso11_comparacion.csv`, `caso11_por_escenario.csv`, `resumen_por_caso.csv`, `comprobacion_reutilizacion.csv`, `sensibilidad_temperatura_plan1.csv`, `plan1_T10_13_por_escenario.csv`. `figuras/`: curvas de Tmax, Tmin S25 y ΔT con límites para M4_SS frío, M4_SS central y M1_SS central (7 casos y 4 temperaturas), vistas 3D con escala común (10–60 °C) a 164 h y 336 h (M4_SS frío) y 210 h (M4_SS central), y la figura de sensibilidad. Los campos 3D completos quedan en `results/bateria7/campos/` (ignorado por git).
