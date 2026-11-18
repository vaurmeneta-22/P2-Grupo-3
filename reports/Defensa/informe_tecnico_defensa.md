# Informe técnico para la defensa (jueves 19-11-2026) — verificación y búsqueda de configuración

> **PROPUESTA FINAL (alternativa B), bajo los supuestos documentados — Anexos 9 y 10.** Colocación continua progresiva a 36 m³/h (D = 17,69 h), 10 °C, manta PUR de 80 mm (4 capas de 20 mm), aislación inferior de 100 mm (propuesta térmica condicionada a viabilidad estructural) y retiro de 20 mm los días 56, 77, 98 y 119 (14-01, 04-02, 25-02, 18-03-2027). Configuración exacta: `configs/defensa/propuesta_B.json`; comparación con la propuesta del compañero (A): `comparacion_propuestas.md`. **Los candidatos de las secciones 1–8 y de los Anexos 1–8 (p. ej. 63,7 m³/h, aislación inferior de 50 mm, retiro cada 14 d) son pasos de la búsqueda y quedan como evidencia; la regla de retiro cada 14 d (Anexo 8, regla A-14d) falla (ΔT 20,71 °C) y se conserva como comparación.** El cumplimiento corresponde a los escenarios y supuestos evaluados (extensión climática supuesta desde el día 28; malla 0,25 m con sensibilidad espacial observada); no es una solución de obra validada.


Elaborado el martes 17-11-2026. No es el artículo ni la presentación. Todo está calculado y disponible sin internet.

## 1. Actualización del repositorio
- Rama `main`. SHA anterior `d9600b45cbc7444e653d579cfe00c9caf14e2ee7`; SHA nuevo `64cef44a6e1fadf6afd73a667d0566e89db5b9c4` (commit de Vicente "Agrega dossier termico final y evidencias"). Fast-forward, sin conflictos (los archivos entrantes no coincidían con ningún archivo local no versionado).
- Respaldo previo: copia local fuera del repositorio (archivos personales sin versionar, estado git y SHA anterior); no se publica. Nada se borró. **Sin commit ni push.**
- Pruebas: 184 pruebas existentes OK (`python -m unittest discover -s tests`, con `PYTHONPATH=src`); con 5 pruebas nuevas (`tests/test_defensa.py`) son 189, todas OK. (pytest no está instalado en `.venv`; se usó unittest.)
- Caso base de Vicente reproducido exactamente (10 h, T colocación 10 °C, PUR 40 mm desde 15 h, suelo 2 m/20 °C, M1_SS, GEFS central, 14 d): Tmax 50,735 °C (208,83 h), ΔT 17,912 °C (336 h). Comando: `python scripts/defensa_busqueda.py configs/defensa/planes_repro.json --curvas M1_SS --climas central`.

## 2. Hallazgo central
**El máximo de ΔT ocurre al final de los 14 días para las configuraciones con aislación, y sigue creciendo después.** Por eso se extendió a 28 y 56 d. Con la manta de 40 mm de Vicente (sin retiro) el ΔT llega a 20,05 °C a los 23,8 d; si la manta se retira a las 375 h (su configuración original: 360 h tras instalar), ΔT = 30,97 °C. Su candidato anterior `propuesto_6h_protegido` pasa a 14 d en M1_SS pero llega a ΔT = 38,8 °C (M1_SS) al retirar la manta a 360 h. Con 14 d solamente, ambos parecen cumplir.

## 3. Criterios y auditoría de S25
- M_T y M_Δ según la definición pedida (`src/mcplan/thermal/decision.py`), registrados en **todos los pasos del solver** (no solo fotogramas), con instante y ubicación. Falla = M_T > 70 °C o M_Δ > 20 °C (una sola falla conjunta si ocurren ambos).
- S25: Δx = 0,25 m no resuelve 25 mm. Implementación auditada: se reconstruye la temperatura de la cara con la condición de borde vigente (Robin con resistencias de moldaje/manta/suelo) y se interpola linealmente hacia el centro del voxel de borde hasta 25 mm. Es una aproximación explícita (error de orden de la curvatura del perfil dentro de media celda). Verificado: la máscara activa se actualiza al activar cada vertida y al instalar/retirar protecciones (excluye aire, suelo y vertidas futuras; la cara inferior en contacto con el suelo cuenta como superficie exterior de hormigón, como en el código existente). Se conservan métricas con nombre propio: `M_delta_centro_voxel` y `M_delta_cara`; además `S25 sin cara inferior` en las corridas de sensibilidad.
- Sensibilidad de S25 a la malla (candidato `f_P2_L2m_U63`, M4_SS fría, ΔT a 56 d): dx 0,5 → 18,35; 0,25 → 17,90; 0,125 → 18,19 °C. Efecto temporal (Δt 600/300/150 s): < 0,01 °C. **Sensibilidad numérica observada ≈ ±0,5 °C (no es una cota rigurosa)** (tabla `sensibilidad_numerica.csv`). Ninguna clasificación del candidato ni de `f_P2_L0_8m_U63` cambia; el competidor más ligero `T10 PUR60 AI50` queda en 19,63 °C (indeterminado), y `T10 PUR80 sin AI` pasa a 20,52 °C (no cumple).

## 4. Datos reales e incertidumbre
**Curvas de hidratación.** Solo existe `data/hidratacion/source_terms.txt` (columnas M1…M4 de dos campañas: SS = 13–20 abril 2025, 22A = 22 abril 2025). **No encontré los archivos originales ni metadatos del profesor** en el repo ni en las carpetas del proyecto. Auditoría (`curvas_auditoria.csv`): ya vienen como Qdot volumétrico en W/m³ (invertidas por el profesor/curso), tiempo en h con origen = inicio de la ventana (reinicia por vertida en el solver), paso 10 min, calor total coincide con el encabezado (±0,06 %), SS dura 240,3 h y 22A 163 h (cola = 0; Q final ≤ 0,3 % del pico). No se deriva ni se multiplica por ρ·cp. Observación: el ensayo usó ρ·cp = 2476·1023 y el modelo 2400·900; con Qdot directo el calor total es el medido, pero el ascenso adiabático del modelo (≈41–46 K SS) difiere del del ensayo (≈35–40 K).
**No se sabe si son réplicas de una misma mezcla.** Podrían ser 4 probetas de una mezcla en 2 campañas o mezclas distintas (el IFC dice M1; 22A es de menor confianza según el archivo). Por eso se calcularon tres interpretaciones de pesos fijadas antes de comparar planes: (a) `nominal_M1_SS` (1 curva); (b) `solo_SS` (4 curvas, 1/4 c/u); (c) `ocho_curvas` (1/8 c/u, enumeración de la pauta, independencia curva-clima). **Pendiente: confirmar con el profesor (o su archivo original) qué representan M1–M4.** Con 24 escenarios en (c) un solo escenario fallido pesa ≈ 4 % (< 5 %), dos ya incumplen.
**Clima.** Sitio oficial lat −31,1115, lon −71,6147. Ensamble NOAA GEFS 0,5° (31 miembros, Open-Meteo, descargado 16-11-2026 21:32 UTC, `data/weather/2026-11-16T213220Z_GEFS_ensemble/original.json`, sha256 en metadatos). t0 común = jueves 19-11-2026 10:00 hora de Ovalle (UTC−3) = 13:00 UTC. Tres trayectorias (frío 11/31, central 10/31, cálido 10/31; medoides de terciles de temperatura media a 14 d, elección de Vicente, conservada sin cambios). El snapshot antiguo que termina el 25-nov no se usa; los 14 d están cubiertos por el nuevo ensamble. **Extensión a 28 d** (nuevo, `scripts/preparar_clima_extendido_gefs.py`): mismos miembros, horas 336–672 también de `original.json` (hasta 17-dic); más allá de ~15 d GEFS no tiene habilidad de pronóstico, es solo una trayectoria del ensamble. **De 28 a 56 d** se repite el último ciclo diario (`extension: ciclo_diario`, etiquetado en el loader). Las probabilidades quedan condicionadas a estas tres historias (peso = fracción de miembros).

## 5. Búsqueda de configuración
Geometría IFC final y volumen fijos (636,8 m³). Supuestos declarados: una cuadrilla, vertidas sin traslape, caudal R; D = fin de la última colocación − t0 (incluye esperas); **activación térmica instantánea al inicio de cada vertida** (el solver la admite así; la duración física sí cuenta en D). Constantes del equipo conservadas: M1 k=2, ρ=2400, cp=900; h=8; moldaje 18 mm k=0,1 en V1/V2/V5 retiro 72 h (IFC); PUR k=0,025 instalado a las 15 h. Controles buscados: temperatura de colocación (10/15/18 °C), espesor de manta PUR (40–100 mm), aislación bajo la base (0/50/100 mm), secuencia/esperas, retiro. Cribado en 128 corridas (`cribado*_resultados_por_escenario.csv`) y comparación completa de 15 planes × 24 escenarios hasta 56 d.

**Paquete que cumple:** `P2` = colocación a 10 °C + manta PUR 80 mm sin retiro + aislación inferior 50 mm (k=0,04). Con **colocación continua** (sin esperas), P_T = 0 y **P_f = 0 en las tres interpretaciones y en 14, 28 y 56 d**, con suelo efectivo de 2 m y de 0,8 m. Máximo ΔT del peor escenario (M4_SS fría): 17,9 °C (suelo 2 m) y 19,1 °C (0,8 m; 19,35 con dx 0,125); máximo ocurre a los ~34–38 d. M_T máx 56,4 °C. Alternativa `P3` (15 °C, PUR 100 mm, aislación 100 mm): ΔT máx 18,7/19,2 °C, M_T 61,3 °C.

Resumen por plan (interpretación `ocho_curvas`; paquete P2 salvo indicación; suelo 2 m salvo indicación; P_f a 14 d / 56 d):

| Plan | D [h] | P_f 14 d | P_f 56 d | máx ΔT 56 d | Cumple (ambos) |
|---|---:|---:|---:|---:|---|
| Continua 100 m³/h (hipotético) | 6,4 | 0 | 0 | 18,4 | Sí |
| **Vertida única continua 63,7 m³/h** | **10,0** | 0 | 0 | 18,4 | **Sí** |
| Vertida única idealizada (todo en t0, D=10 h) | 10,0 | 0 | 0 | 18,3 | Sí (cota) |
| Continua 36 m³/h (una cuadrilla) | 17,7 | 0 | 0 | 18,4 | Sí |
| Una cada 4,5 h (más rápida que propuesta) | 22,5 | 0 | 0 | 18,4 | Sí |
| Propuesto 6 h (IFC) | 30 | 0,125 | 0,125 | 20,7 | **No** (1 curva) |
| Cada 8 h (más lento) | 40 | 0,875 | 0,875 | 24,5 | No |
| Cada 12 h (más lento) | 54 | 1,0 | 1,0 | 28,5 | No |
| Continua 63,7, suelo 0,8 m | 10,0 | 0 | 0 | 19,4 | Sí (margen 0,6) |
| Continua 63,7, paquete P1 (10 °C, PUR 40, sin AI) | 10,0 | 0,169 | 0,415 | 23,5 | No (nominal M1_SS: 20,05, indeterminado) |
| `ref_protegido_6h` (candidato anterior) | 30 | 0,375 | 1,0 | 43,8 | No |

Los planes escalonados fallan porque, por definición, S25 incluye el hormigón recién colocado (10 °C) frente al núcleo caliente de las capas previas: esperas largas penalizan. En el plan de 6 h el resultado depende de cómo se representa la colocación (ΔT del caso crítico M4_22A cálida: 19,42 °C instantánea, 21,46 °C progresiva, 19,24 °C al final), por lo que está en la frontera.

**Estado de factibilidad / recomendación.** El menor D factible entre los planes evaluados es la colocación continua; con 63,7 m³/h (caudal implícito de la secuencia de 10 h de Vicente) D = 10 h, con 36 m³/h D = 17,7 h, y 6,4 h con 100 m³/h (hipotético). **D queda fijado por el caudal disponible, no por la térmica**: la cantidad de bombas/plantas que sostenga 63,7 m³/h (V1 y V2 exigen ≈ 74 m³/h instantáneos) no está validada. Recomendación: **vertida continua + paquete P2**. No se afirma optimalidad global (búsqueda acotada a la familia descrita) ni ausencia de fisuración.

## 6. Sensibilidades de supuestos (`sensibilidad_fisica.csv`, candidato, ΔT a 56 d; M4_SS fría / M1_SS fría; referencia 17,90 / 15,19)
- **Retiro de manta: decisivo.** Retiro a 14 d (351 h): 36,4/31,4 °C; a 28 d: 32,2/27,6 °C (también incumple). **El día 28 NO es una fecha segura de retiro.** El cumplimiento encontrado corresponde únicamente a mantener la manta durante todo el horizonte simulado de 56 d; no hay ninguna regla de retiro simulada que cumpla (pendiente, ver comprobaciones posteriores).
- k de la manta 0,035 (vs 0,025): 19,67/16,73 (aún cumple, margen estrecho); 0,022: 17,36.
- Suelo: L=0,8 m: 19,06/16,18; L=4 m: 17,23; T suelo 16 °C: 18,66; 24 °C: 17,25. **El suelo efectivo de 2 m no tiene fundamento documentado en lo que revisé (E10 usa 0,8 m); se deja como sensibilidad, y el candidato cumple con ambos.**
- h aire 5/12 W/m²K: ±0,1. Colocación 8/12 °C: 16,93/18,86 (cada °C ≈ 0,24 °C de ΔT). Hora de inicio (+6, +12 h): −0,06. Inicio de manta a 6 o 24 h: ±0,08. Moldaje retirado a 168 h o sin moldaje: ~0.
- Los parámetros h=8, moldaje 18 mm k=0,1, PUR k=0,025 siguen siendo supuestos del equipo no verificados con un producto específico; su efecto es pequeño salvo k de la manta y el retiro.
- Supuestos que **no** se resolvieron: temperatura de colocación 10 °C exige hormigón enfriado (método y costo no justificados en el repo); aislación bajo la base requiere verificación estructural; no se modeló el relleno/entorno real ni la maduración dependiente de la temperatura de la curva.

## 7. Verificación
- Balance de energía (`balance_energia.csv`, caso candidato M4_SS fría a 28 d): calor liberado del solver vs Σ V·∫Q dt esperado: +0,005 %; cierre E_final = E_act + E_fuente + E_frontera con residuo relativo 6×10⁻¹⁶; 23,5 % del calor generado salió al entorno. Las 490 corridas del caché tienen residuo relativo < 2×10⁻¹⁵.
- Convergencia: ver §3. Pruebas: 189 OK.

## 8. Archivos y reproducción (todos locales)
Tablas: `reports/Defensa/tablas/` — `comparacion_planes.csv` (por plan e interpretación: D, P_T, P_Δ, P_f, cumplimiento a 14/28/56 d), `resultados_por_escenario.csv` (una fila por plan×curva×clima, pesos, máximos, instantes y ubicaciones, a 14/28/56 d), sensibilidades, balance, auditoría de curvas. Figuras: `reports/Defensa/figuras/` (curvas de hidratación y clima, cronogramas, series con líneas 70 y 20 °C, D vs P_f con línea 5 %, cortes 3D con escala común 10–60 °C). Corridas crudas en `results/defensa/cache` (ignoradas por git).
```
$env:PYTHONPATH="src"; $env:PYTHONIOENCODING="utf-8"
python scripts/preparar_clima_extendido_gefs.py
python scripts/auditar_curvas_hidratacion.py
python scripts/defensa_generar_planes.py final
# un caso (candidato, peor escenario):
python scripts/defensa_busqueda.py configs/defensa/planes_final.json --solo f_P2_L2m_U63 --curvas M4_SS --climas fria --horizonte-h 1344
# comparación completa (≈70 min con 15 procesos; reanuda sin repetir lo ya calculado):
python scripts/defensa_busqueda.py configs/defensa/planes_final.json --curvas todas --climas todas --procs 15 --horizonte-h 1344
python scripts/defensa_busqueda.py configs/defensa/planes_final2.json --curvas todas --climas todas --procs 15 --horizonte-h 1344
python scripts/defensa_consolidar.py --prefijo f_ --salida results/defensa/tabla_final
python scripts/defensa_figuras.py pareto results/defensa/tabla_final/comparacion_planes.csv
python scripts/defensa_sensibilidad.py ; python scripts/defensa_balance_energia.py f_P2_L2m_U63 M4_SS fria 672
```
Nota: la clave de reanudación incluye el hash del código; si se edita `defensa_busqueda.py` o el solver, se recalcula.

## 9. Pendientes críticos para el jueves
1. Confirmar qué son las 8 curvas (réplicas vs mezclas) y obtener los archivos originales/metadatos del profesor.
2. Justificar o cambiar: colocación a 10 °C (enfriamiento), aislación bajo la base, suelo efectivo (0,8 m vs 2 m, T=20 °C) y caudal de 63,7 m³/h.
3. Definir y simular una regla de retiro gradual de la manta (cualquier retiro simulado, incluso a 28 d, rompe el criterio; ninguna regla cumple aún comprobada) y el cierre real tras 56 d (relleno).
4. Decidir si el criterio oficial es 14 d (la mayoría de paquetes cumple) o el horizonte extendido (solo P2/P3 cumplen); el informe reporta ambos.
5. Refrescar el pronóstico GEFS cerca del jueves si se quiere (`descargar_snapshot_gefs_ensemble.py`); la matriz actual es del 16-11.
6. Decidir qué commitear (nada se commiteó): `scripts/defensa_*`, `tests/test_defensa.py`, `configs/defensa/`, `reports/Defensa/`, cambios en `src/mcplan/thermal/decision.py` y `data/weather/…_28d/`.

---
# Anexo — comprobaciones posteriores (17-11-2026, tarde)
31 corridas, sin errores, ~17 min de reloj con 12 procesos (cada corrida: 2,3–2,6 min las de 56 d, 5,4–7,5 min las de 150 d). Resultados en `results/defensa/cache2/` (reanudable: `python scripts/defensa_comprobaciones.py --procs 12 [prog|retiro|sinAI]`) y tablas `tablas/comprob_progresiva.csv`, `comprob_retiro.csv`, `comprob_sin_AI.csv`. **Solo 3 escenarios críticos (M4_SS fría, M4_SS cálida, M1_SS central): sirven para mecanismos y selección de candidatos, no sustituyen la evaluación completa de P_f.**

**Extensión climática de las corridas a 150 d (supuesto de exploración).** Horas 0–672: trayectorias GEFS reales del miembro (descarga 16-11). Horas 672–3600: continuación por repetición del último ciclo diario de 24 h del pronóstico (método `ciclo_diario` de `mcplan.weather.snapshots.extender_serie`; cada punto queda etiquetado `extension:ciclo_diario`). No es pronóstico.

**1. Colocación progresiva (caudal constante, llenado por capas z↑, instantes cuantizados a 5 min) vs instantánea.** ΔT a 56 d cambia ≤ 0,006 °C en vertida continua a 63,7 y 36 m³/h; ΔT a 14 d: ≤ 0,003 °C (63,7) y +0,25–0,28 °C (36); M_T +0,04–0,05 °C (63,7) y +0,08–0,10 °C (36); D no cambia (10,0 h y 17,7 h; es el fin físico de colocación). Peor caso progresivo: 18,42 °C (M4_SS cálida, 36 m³/h). Caudal: por construcción, 63,7 y 36,0 m³/h con volúmenes IFC por vertida. **Salvedad:** el modelo voxelizado tiene volúmenes por vertida distintos al IFC (V1/V2 −10 %, V3 +16 %; total −0,5 %), por lo que el caudal efectivo del modelo oscila 55–77 m³/h (nominal 63,7) y 31–42 m³/h (nominal 36). Frente a la capacidad declarada por el equipo (36 m³/h, `evaluacion_final.jsonc`), **el plan de 10 h la excede (63,7)**; el de 17,7 h la iguala. La capacidad realista sigue sin confirmar.

**2. Reglas de retiro (4 capas de 20 mm de PUR retiradas de a una en todas las vertidas, 150 d; manta de 80 mm).** Fechas desde t0 de las 4 retiradas; "máx. post-retiro" = máximo de ΔT desde el primer retiro:
| Regla | Retiros [d] | M4_SS fría | M4_SS cálida | M1_SS central | ¿Cumple en los 3? |
|---|---|---|---|---|---|
| sin retiro (base) | — | 17,9 | 18,4 | 15,3 | sí |
| R84_c14d | 84, 98, 112, 126 | 14,3 (a 84,7 d) | 15,1 (a 84,9 d) | 12,3 | sí; máximo global sigue siendo el de ~38–41 d |
| R56_c14d | 56, 70, 84, 98 | 17,2 (57,8 d) | **18,67 (99,8 d)** | 15,1 | sí, margen ≥1,3 °C |
| R56_c7d | 56, 63, 70, 77 | 21,2 | 23,4 | 19,0 | **no** |
| R42_c7d | 42, 49, 56, 63 | 24,1 | 26,4 | 21,4 | **no** |
Candidatas ejecutables: R84_c14d y R56_c14d (parámetros explícitos: 4 capas de 20 mm, intervalo 14 d, inicio 56 u 84 d). Alternativa basada en temperatura (retirar una capa solo si ΔT medido < 15 °C) **no se simuló**. Al final de 150 d el núcleo aún está a 22–29 °C (aire ~18 °C). No incluye relleno ni lluvia/viento. Retiros más rápidos que 14 d entre capas o antes del día 56 incumplen.

**3. Sin aislación inferior (peor escenario M4_SS fría, 56 d; ΔT):**
| Paquete | suelo 2 m | suelo 0,8 m |
|---|---|---|
| T10, PUR 100 / 120 / 150 mm | 19,68 / 19,34 / 19,06 | 24,70 / 24,58 / 24,56 |
| T15, PUR 120 / 150 mm | 21,94 / 21,62 | 27,85 / 27,83 |
Prometedor solo `T10 + PUR 120–150 mm`, **y solo si el suelo efectivo es 2 m**; con 0,8 m ninguna espesor ayuda (la base domina). Márgenes de 0,7–0,9 °C ≈ la sensibilidad numérica observada (indeterminado hasta refinar malla). A 15 °C no hay solución sin aislación inferior en lo probado.

**4. Candidatos que requieren evaluación completa (24 escenarios, ambos suelos, y dx = 0,125) antes de recomendarse:** (a) vertida continua + P2 + retiro R84_c14d o R56_c14d (también 150 d); (b) `T10 + PUR 150 sin aislación inferior`, suelo 2 m; (c) el plan de 17,7 h con colocación progresiva, por ser el que respeta 36 m³/h. P2 con manta fija sigue siendo el único con evaluación completa (ver §5).

**5. Comparación 10 °C vs 15 °C (resultados ya evaluados, 24 escenarios, máx ΔT a 56 d):** 
| Paquete | Medidas exigidas | suelo 2 m | suelo 0,8 m | M_T máx |
|---|---|---|---|---|
| P2: 10 °C | enfriar el hormigón a 10 °C; PUR 80 mm; aislación inferior 50 mm | 18,4 (margen 1,6) | 19,4 (0,6) | 56,4 |
| P3: 15 °C | PUR 100 mm; aislación inferior 100 mm | 18,7 (1,3) | 19,2 (0,8) | 61,3 |
Ambos exigen aislación inferior; 15 °C reduce la exigencia de enfriamiento pero duplica el espesor bajo la base. Ninguno se eligió por aprobación: son los dos paquetes del cribado que superan los 24 escenarios.

**6. Datos a confirmar**
| Dato | Qué falta | Quién | Efecto en la decisión |
|---|---|---|---|
| M1–M4 y campañas SS/22A | ¿réplicas de una mezcla o mezclas distintas? Archivos originales/metadatos del ensayo, mezcla del IFC (M1) | Profesor / Vicente | Define pesos y P_f (3 interpretaciones hoy) |
| Capacidad de colocación | m³/h sostenibles, nº de bombas; ¿36 (declarado) o 63,7? | Equipo / constructora | Define D (10 vs 17,7 h) |
| Temperatura de colocación | ¿10 °C alcanzable (hielo/agua fría/áridos)? ¿15 °C? | Equipo / proveedor | Decide P2 vs P3 |
| Aislación inferior | Viabilidad estructural y constructiva (capacidad de apoyo), producto y k | Equipo / ingeniería | Sin ella solo 'T10+PUR≥120, suelo 2 m' (marginal) |
| Manta | Producto, k real (0,025 supuesto; 0,035 → +1,8 °C), forma de retiro por capas | Equipo / proveedor | Cumplimiento y regla de retiro |
| Suelo efectivo | 0,8 m vs 2 m y T del suelo | Equipo | Margen 1,6 vs 0,6 °C |

Las conclusiones siguen condicionadas a: interpretación de curvas y pesos, 3 historias GEFS (frío/central/cálido, medoides del 16-11), extensión climática a 28 y 56+ d, y activación instantánea (cotejada con progresiva). La sensibilidad numérica ±0,5 °C es observada (dx 0,5/0,25/0,125), no una cota.

---
# Anexo 2 — cierre centrado en el candidato (18-11-2026)
**Candidato:** colocación continua **progresiva** a 36 m³/h (una cuadrilla; D = 17,7 h), 10 °C, PUR 80 mm (desde las 15 h, sin retiro dentro de 56 d), aislación inferior 50 mm. 48 corridas (24 escenarios × suelos 2 m y 0,8 m; 14/28/56 d); ~3 min cada una; sin errores. Tablas: `comprob_candidato_planes.csv`, `comprob_candidato_escenarios.csv`.

| Suelo | Interpretación de pesos | P_T | P_Δ | P_f (14/28/56 d) | Máx ΔT 14 d / 28 d / 56 d [°C] | Máx M_T [°C] |
|---|---|---|---|---|---|---|
| 2 m | M1_SS · 4 SS · 8 curvas | 0 | 0 | 0 / 0 / 0 en las tres | 16,4 / 17,6 / 18,4 (8 curvas); M1_SS: 13,1 / 15,0 / 15,7 | 56,6 |
| 0,8 m | M1_SS · 4 SS · 8 curvas | 0 | 0 | 0 / 0 / 0 en las tres | 16,4 / 19,1 / 19,4 (8 curvas); M1_SS: 13,4 / 16,2 / 16,6 | 56,5 |
Peor escenario: M4_SS cálida (ΔT 18,42 °C a 41,8 d con suelo 2 m; 19,45 °C a 36,9 d con 0,8 m), en la base (x≈−9, y≈−3,6, z≈0,02–0,12 m). **Suelo 0,8 m: cumple con margen de 0,55 °C**; la sensibilidad numérica observada en un caso análogo (dx 0,125 sumó +0,29 °C) reduciría ese margen a ≈0,3 °C: *cumple, pero sin robustez*; no se corrió dx = 0,125 para este candidato. Los máximos de ΔT ocurren después del día 14 (respaldado por pronóstico solo hasta ~día 14; la ventana 14–28 d usa trayectorias GEFS sin habilidad; de 28 a 56 d el ciclo diario repetido). **Conclusiones respaldadas por datos meteorológicos: 14 d; las de 28 y 56 d dependen de la extensión climática declarada.**

**Volúmenes por vertida (malla vs IFC)** — `volumenes_ifc_vs_voxel.csv`. A dx = 0,25: V1 y V2 −9,8 %, V3 +16,0 %, V4 −10,0 %, V5 +0,3 %, total −0,53 %; con otras mallas cambian de signo (dx 0,5: V4 +37,8 %; dx 0,125: V2 +13,3 %): son artefactos de a qué vertida se asigna cada voxel en las interfaces (alturas 0,55/1,10/2,20 m no alineadas con la malla), no convergen al refinar. Hay que distinguir: (a) **caudal constructivo según IFC** = V_IFC / duración = 36,0 m³/h constante, D = 17,69 h; (b) **volumen activado numéricamente**, que a dx = 0,25 produce 31–42 m³/h por hora con el mismo cronograma. Para medir el efecto se repitió la progresiva con duraciones = volumen numérico / 36 (caudal numérico constante, D = 17,59 h): ΔT cambia ≤ 0,009 °C y M_T ≤ 0,02 °C en 4 escenarios (`vol_*` en `cache2`). Efecto residual de la asignación de volumen sobre el calor liberado por vertida y la geometría: no aislado; su cota observada es la dispersión entre mallas (≈0,45 °C en ΔT). Instantánea vs progresiva en los 24 escenarios con suelo 2 m: ΔT a 56 d +0,003 a +0,23 °C, a 14 d ≤ +0,29 °C, M_T ≤ +0,11 °C.

**Opción de 10 h (63,7 m³/h): solo sensibilidad de mayor capacidad**, no recomendación (excede los 36 m³/h declarados). Resultados existentes: P_f = 0 en las tres interpretaciones, ΔT máx 18,4 (suelo 2 m) / 19,4 (0,8 m).

**Retiro de manta (20 mm cada 14 d desde el día 56) — pendiente para el candidato completo.** Ya probado solo en 3 escenarios con vertida continua a 63,7 m³/h (máx post-retiro 18,67 °C a los 99,8 d). Costo previsto para comprobar el candidato completo hasta 150 d (3600 h): ~6,5 min por corrida con 12 procesos; 24 escenarios con suelo 2 m ≈ 13–15 min; con ambos suelos (48) ≈ 26–30 min; RAM moderada. Conclusiones de esas corridas dependen de la extensión climática (ciclo diario repetido desde la hora 672).

**Cambios en la recomendación:** pasa a ser el candidato de 17,7 h progresivo (recursos del equipo); el de 10 h queda como sensibilidad. Se mantienen todas las comparaciones anteriores.

**Decisiones pendientes con el equipo/profesor**
| Tema | Pregunta concreta | Por qué importa |
|---|---|---|
| Ocho curvas | ¿M1–M4 y SS/22A son réplicas de una mezcla o mezclas distintas? ¿Cuál es la mezcla del IFC? Archivos originales/metadatos | Define pesos y P_f (hoy 3 interpretaciones, todas dan P_f = 0 para el candidato) |
| Colocar a 10 °C | ¿Se puede enfriar el hormigón a 10 °C con el aire a 12–25 °C? ¿Método, costo, control en obra? | Cada °C ≈ 0,24 °C de ΔT; a 15 °C se exige más aislación (PUR 100 + 100 mm bajo la base) |
| Aislación inferior | ¿Es viable bajo la base (resistencia a compresión, apoyo, drenaje)? Producto y k (0,04 supuesto) | Sin ella el candidato no cumple (≥ 19,1–24,7 °C según el suelo) |
| Manta | Producto, espesor real, k (0,025 supuesto; con 0,035 el ΔT sube ~1,8 °C), forma de retirarla por capas | Margen de ΔT y regla de retiro |
| Suelo | ¿Qué fundamenta 2 m, 0,8 m y T = 20 °C? | Margen 1,6 °C (2 m) vs 0,55 °C (0,8 m) |

---
# Anexo 3 — retiro de manta del candidato completo y volúmenes por vertida (18-11-2026)
## A. Retiro: una capa de 20 mm cada 14 d desde el día 56 (retiros a 56, 70, 84, 98 d), candidato 36 m³/h progresivo
48 corridas (24 escenarios × 2 suelos) hasta 150 d (3600 h), ~7 min cada una con 12 procesos, ~32 min en total, sin errores. Tablas: `comprob_retiro_candidato_escenarios.csv` (por escenario: máximos de T y ΔT a 14 d, 56 d y 150 d con instante y ubicación; máximos tras cada retiro y tras el retiro completo) y `comprob_retiro_candidato_planes.csv`.
**Dependencia climática:** las horas 0–336 son pronóstico GEFS; 336–672 son trayectorias GEFS sin habilidad; **672–3600 h repiten el último ciclo diario (supuesto de exploración, `extension:ciclo_diario`)**. Todo lo posterior al día 28, y en particular los retiros (días 56–98) y sus picos, depende de ese supuesto; las conclusiones respaldadas por pronóstico siguen siendo las de 14 d.
| Suelo | P_T | P_Δ | P_f (56 d → 150 d) en las 3 interpretaciones | Máx ΔT 56 d → 150 d (8 curvas) | Máx M_T |
|---|---|---|---|---|---|
| 2 m | 0 | 0 | 0 → 0 | 18,42 → **18,73 °C** | 56,6 °C |
| 0,8 m | 0 | 0 | 0 → 0 | 19,45 → 19,45 °C | 56,5 °C |
Máximo de ΔT posterior a cada retiro (peor escenario de los 24, valor / suelo 2 m · suelo 0,8 m): tras retiro 1 (día 56): 18,10 · 18,25; tras 2 (día 70): 17,86 · 17,24; tras 3 (día 84): 17,55 · 16,26; tras 4 (día 98, retiro completo): **18,73 · 16,77 °C**. Con suelo 2 m el máximo global a 150 d (M4_SS cálida, 18,73 °C a los 99,8 d, cara lateral x≈+9, y≈+3,6, z≈1,1 m) ocurre ~2,4 días después del último retiro y supera en 0,31 °C el máximo previo (18,42 °C); su margen es 1,27 °C. Con suelo 0,8 m el máximo (19,45 °C, día 36,9, base) es anterior a los retiros: los retiros no lo empeoran. Excesos del máximo a 150 d sobre el de 56 d: ≤0,73 °C (suelo 2 m, 21 escenarios) y ≤0,01 °C (0,8 m). Máximos por curva a 150 d (2 m): M4_SS 18,73; M2_SS 17,88; M3_SS 17,58; M1_SS 16,46; 22A ≤16,4. Ningún escenario falla. **Tendencias al final (días 140–150):** ΔT decrece (−0,05 a −0,14 °C/día), ΔT final ≤ 7,3 °C; Tmax del núcleo 21–24 °C y bajando hacia el aire (~18 °C). No hay crecimiento que obligue a extender más el horizonte para este candidato bajo la hipótesis climática. La ubicación de los máximos posteriores a cada retiro no se almacenó (solo la del máximo global).

## B. Volúmenes por vertida vs dx (corrige la afirmación anterior «no convergen»)
Tabla `volumenes_convergencia_dx.csv` (error absoluto en m³, voxel − IFC):
| malla | V1 | V2 | V3 | V4 | V5 | total |
|---|---:|---:|---:|---:|---:|---:|
| 0,5 | −13,1 | −13,1 | −2,2 | +33,6 | +0,6 | +5,7 |
| 0,4 | −41,3 | +65,0 | −52,8 | +25,2 | +7,9 | +3,9 |
| 0,3 | +13,7 | +13,7 | −40,8 | +11,8 | −3,4 | −4,9 |
| **0,25** | −14,5 | −14,5 | +34,3 | −8,9 | +0,1 | −3,4 |
| 0,2 | ERROR de voxelización (solape V4/V5) |
| 0,15 | +12,8 | −27,3 | +21,1 | −6,2 | +1,7 | +2,0 |
| 0,125 | −13,8 | +19,6 | −0,3 | −8,7 | 0,0 | −3,2 |
| 0,1 | ERROR (solape V1/V2) |
| 0,075 | −6,9 | +13,2 | −10,3 | +3,1 | −1,0 | −1,8 |
| 0,0625 | +3,3 | +3,3 | −8,3 | +2,1 | 0,0 | +0,3 |
| 0,25×0,25×**0,05** | −1,2 | −1,2 | +0,04 | +0,09 | +0,12 | −2,1 |
| 0,125×0,125×0,05 | −0,5 | −0,5 | −0,5 | −1,1 | 0,0 | −2,5 |
**Diagnóstico.** (1) La afirmación anterior era imprecisa: los errores absolutos **sí disminuyen con dx, pero de forma no monótona y lenta** (la magnitud máxima por vertida pasa de ~65 m³ a ~8 m³ entre 0,4 y 0,0625 m). (2) La causa principal es **geométrica**: los límites entre vertidas están en z = 0,55; 1,10; 2,20; 3,30 y 4,30 m y, con muestreo por centro de voxel, cada interfaz se redondea al plano de malla más cercano; el error de volumen ≈ (desplazamiento de la interfaz) × (área de la sección) (p. ej. a dx = 0,25: V1 queda con 0,50 m en vez de 0,55 m → −14 m³; V3 con 1,25 m en vez de 1,10 m → +34 m³). Con una malla alineada en z (dz = 0,05 m) los errores por vertida caen a ≤ 1,2 m³ (≤ 0,8 %); el resto (−2 m³) es la escalera lateral del contorno octogonal/circular. (3) Hay además un **problema de asignación en las interfaces**: el IFC deja un solape de ~3–6×10⁻⁸ m entre vertidas (p. ej. V4 termina en 3,30000001 y V5 empieza en 3,29999995), y el voxelizador aborta cuando un centro de voxel cae en esa franja (dx = 0,2 con z = 3,3 exacto; dx = 0,1 y dz = 0,1 con z = 0,55). No afecta a las mallas usadas (0,5; 0,25; 0,125) pero impide barridos de dx más finos con ciertos pasos; corrección propuesta (no aplicada): regla de desempate determinista en `voxelizar_vertidas` (asignar al centro a la vertida inferior/anterior).
**Corrección geométrica que se propone ensayar:** malla anisótropa 0,25 × 0,25 × 0,05 m (interfaces alineadas), que cuesta ≈ 5× celdas.

## C. Comprobación espacial del candidato exacto con suelo 0,8 m (propuesta de costo; no lanzada)
Escenario: M4_SS cálida (máximo ΔT = 19,45 °C), progresivo 36 m³/h, retiro no necesario (el máximo es anterior), horizonte 56 d (1344 h). Opciones, una corrida cada una, ejecutadas solas (la concurrencia no aporta):
1. Isótropa dx = 0,125 (≈ 8× celdas): ≈ 10–15 min, ~3 GB de RAM (en el lote anterior, 8 en paralelo tardaron ~39 min cada una).
2. Anisótropa 0,25×0,25×0,05 (interfaces alineadas; S25 a 25 mm coincide con el centro del voxel de borde en z): ≈ 5–10 min, ~1,5 GB (paso temporal limitado por dz ≈ 360 s).
3. Combinada 0,125×0,125×0,05 (~19× celdas): ≈ 30–45 min, ~5–6 GB.
Recomendación: 1 y 2 (≈ 20–25 min en total); 3 solo si 1 y 2 discrepan en más de 0,3 °C. El margen de 0,55 °C se considera comparable con la sensibilidad entre mallas hasta tenerlas.

---
# Anexo 4 — comprobación espacial del candidato exacto, suelo 0,8 m, M4_SS cálida (18-11-2026)
Candidato: progresivo 36 m³/h, 10 °C, PUR 80 mm (desde 15 h, sin retiro), aislación inferior 50 mm, horizonte 56 d, dt máx. 600 s. Script `scripts/defensa_malla_candidato.py`; resultados `results/defensa/cache2/malla_*.json` y `reports/Defensa/tablas/comprob_malla_candidato.json`.
**Verificación previa de la malla anisótropa 0,25×0,25×0,05 (73×73×86):** el solver admite `dx_m` distinto por eje (ya cubierto por `test_extremos_termicos`, dx_m=(0,5; 0,4; 0,3)); volumen de celda 0,003125 m³; áreas de cara 0,0125 (x, y) y 0,0625 m² (z); paso estable por conductancias 1250 s > 600 s, así que el paso real es el mismo (600 s) que en la referencia (en dx = 0,125 el estable es 2812 s); S25 usa d por eje: fracción 0,2 en x/y y 1,0 en z (el punto a 25 mm coincide con el centro del voxel de borde, profundidad = d/2, admitida por la guardia); balance de energía cerrado en las tres corridas (residuo ≤ 8×10⁻⁶ J); energía generada del solver vs Σ volumen·∫Q: +0,067 % (dx 0,125) y +0,066 % (aniso).
| Malla | ΔT máx 56 d | instante / ubicación S25 mín | Tmax (instante) | margen a 20 °C | ΔT 14 d | cumple | tiempo |
|---|---:|---|---|---:|---:|---|---:|
| ref 0,25³ | 19,446 | 885,7 h; (−8,88; −3,63; 0,025) | 56,50 (223,3 h) | 0,554 | 15,39 | sí | — |
| 0,125³ | 19,740 (+0,29) | 860,7 h; (−8,94; −3,69; 0,025) | 56,57 (222,8 h) | 0,260 | 15,74 | sí | 10,7 min |
| 0,25×0,25×0,05 | 19,842 (+0,40) | 882,7 h; (−8,98; −3,63; 0,025) | 56,64 (222,3 h) | **0,158** | 15,80 | sí | 6,3 min |
Volumen y energía generada por vertida (V1…V5; vol m³ voxel, error vs IFC, energía GJ): ref (0,25³): 133,1/133,1/248,6/80,1/38,5 (−14,5/−14,5/+34,3/−8,9/+0,1); 0,125³: 133,8/167,2/213,9/80,3/38,4 (−13,8/+19,6/−0,3/−8,7/0,0), energía 13,4/16,8/21,4/8,0/3,8; aniso: 146,4/146,4/214,3/89,0/38,5 (−1,2/−1,2/+0,04/+0,09/+0,12), energía 14,67/14,67/21,46/8,92/3,86. La malla anisótropa reproduce los volúmenes IFC por vertida (≤0,8 %); las otras no.
**Lectura.** Las dos refinamientos aumentan ΔT (no lo reducen), la corrida ref estaba subestimada y la malla con mejor geometría da 19,84 °C: **el margen cae a 0,16 °C, menor que la diferencia entre mallas (0,29–0,40 °C)**. No hay convergencia demostrada (0,25→0,125: +0,29; 0,25→dz 0,05: +0,40) y si ambos efectos (resolución horizontal y alineación vertical) se sumaran, la combinada podría superar 20 °C. **Con suelo 0,8 m el cumplimiento queda INDETERMINADO**, aunque las tres mallas probadas den < 20 °C. Con suelo 2 m el máximo a 150 d es 18,73 °C (margen 1,27 °C): sumar +0,3–0,4 °C lo mantendría ≥ 0,8 °C bajo el límite (extrapolación, no comprobada con malla fina). La tercera corrida (0,125×0,125×0,05, ~30–45 min, 5–6 GB) **no se lanzó**.
**Alcances:** P_f = 0 corresponde al conjunto discreto evaluado (3 trayectorias GEFS × 8 curvas con las interpretaciones de pesos declaradas), a malla 0,25 m; el retiro está comprobado solo bajo la extensión climática supuesta (ciclo diario repetido desde la hora 672); el desempate del voxelizador no se aplicó. No hay commit ni push.

---
# Anexo 5 — tercera comprobación espacial (malla 0,125×0,125×0,05 m), candidato exacto, suelo 0,8 m, M4_SS cálida, 56 d
Misma configuración y mismos parámetros físicos; dt máximo 600 s en las cuatro mallas (paso estable 1250 s en la anisótropa, 1023 s en la combinada, 2812 s en 0,125³: nunca limita), por lo que la comparación es espacial y no temporal (sensibilidad temporal previa < 0,01 °C). Corrida sola: 27,8 min, ~19 GB de RAM (campo guardado cada hora; mayor que los 5–6 GB estimados).
| Malla | voxeles | ΔT máx 56 d | instante / ubicación S25 mín (x,y,z) | Tmax (instante; ubic.) | margen a 20 °C | balance [J] / fuente vs esperada |
|---|---:|---:|---|---|---:|---|
| 0,25³ | 40 536 | 19,446 | 885,7 h; (−8,88; −3,63; 0,025) | 56,50 (223,3 h; −0,13, 0,13, 2,38) | 0,554 | — |
| 0,125³ | 324 380 | 19,740 | 860,7 h; (−8,94; −3,69; 0,025) | 56,57 (222,8 h; −0,06, −0,06, 2,31) | 0,260 | −7,6×10⁻⁶ / +0,067 % |
| 0,25×0,25×0,05 | 203 088 | 19,842 | 882,7 h; (−8,98; −3,63; 0,025) | 56,64 (222,3 h; −0,13, −0,13, 2,23) | 0,158 | −3,8×10⁻⁶ / +0,066 % |
| **0,125×0,125×0,05** | 811 914 | **19,957** | 859,8 h; (−8,98; −3,69; 0,025) | 56,64 (222,3 h; −0,06, −0,06, 2,23) | **0,043** | −3,8×10⁻⁶ / +0,066 % |
Volúmenes por vertida (m³; error vs IFC) en la combinada: 147,16 (−0,47) / 147,16 (−0,47) / 213,77 (−0,46) / 87,85 (−1,08) / 38,38 (0,00); energía generada 14,74/14,74/21,41/8,80/3,84 GJ.
**Conclusión.** El ΔT crece de forma monótona con el refinamiento (19,45 → 19,74 → 19,84 → 19,96 °C); los efectos horizontal y vertical se suman casi linealmente. En la malla más fina el máximo queda **por debajo de 20 °C por solo 0,043 °C**: *no* incumple en esta comprobación, pero **no queda demostrado un cumplimiento robusto**, y la tendencia monótona hace plausible que mallas aún más finas lo superen (no comprobado). Para suelo 0,8 m el candidato debe tratarse como **no confirmado / sin margen**. Esta corrida verifica un único escenario crítico a un solo suelo; no reemplaza la evaluación completa de probabilidades en esta malla. Con suelo 2 m el margen de 1,27 °C (máx 18,73 °C a 150 d) se reduciría en ≈ 0,5 °C por este sesgo de malla (extrapolación) y seguiría ≳ 0,7 °C.
**Ajuste físico acotado propuesto (no lanzado):** aumentar la aislación inferior de 50 a 100 mm (k = 0,04), sin tocar el suelo. En el cribado (malla 0,25³, M4_SS fría, 56 d, 10 °C, PUR 80 mm) pasar de AI 50 a 100 mm bajó ΔT de 18,89 a 16,94 °C con suelo 0,8 m (−1,95 °C) y de 17,90 a 17,24 °C con suelo 2 m. Costo de comprobarlo: (a) escenario crítico en malla 0,25³: ~3 min; (b) en 0,25×0,25×0,05: ~6 min; (c) en la malla fina 0,125×0,125×0,05: ~28 min y ~19 GB; (d) si se acepta, evaluación completa de 24 escenarios × 2 suelos en malla 0,25³ (progresivo, 56 d): 48 corridas, ~15–18 min con 12 procesos, y retiro hasta 150 d: ~32 min. Recomendación de orden: (a)+(b)+(c) ≈ 40 min, luego (d). Requiere además verificar la viabilidad estructural de 100 mm bajo la base.

---
# Anexo 6 — aislación inferior de 100 mm (propuesta térmica condicionada a viabilidad estructural y constructiva)
Candidato idéntico (progresivo 36 m³/h, 10 °C, PUR 80 mm, suelo 0,8 m sin cambios) con aislación inferior de **100 mm** (k = 0,04 supuesto) en vez de 50 mm. Escenario M4_SS cálida, 56 d, dt máx 600 s. Script: `python scripts/defensa_malla_candidato.py correr <ref|aniso005|comb> 100 reducido`; resultados `results/defensa/cache2/malla_*_AI100_red_*.json`.
**Almacenamiento reducido** (`correr_reducido`): los observadores calculan M_T y M_Δ en cada paso; solo se guardan campos completos en t = 0, 336 h, 860 h y final. Validación: candidato AI 50, malla 0,25³: ΔT 19,445878 °C (reducido) vs 19,445865 °C (muestreo horario), Tmax 56,49731 vs 56,49729 °C (diferencia ≤ 1,3×10⁻⁵ °C). RAM de la malla fina: ~4 GB (antes ~19 GB). Ejecución sin concurrencia.
| Malla | AI 50 mm: ΔT máx | AI 100 mm: ΔT máx | instante / ubicación S25 mín (x,y,z) | Tmax AI 100 (instante) | margen a 20 °C | balance [J] | tiempo |
|---|---:|---:|---|---|---:|---:|---:|
| 0,25³ | 19,446 | **18,157** | 1051,5 h; (−8,98; −3,63; 1,125) lateral | 56,56 (224,1 h) | 1,84 | −5,3×10⁻⁵ | 0,5 min |
| 0,25×0,25×0,05 | 19,842 | **18,460** | 1027,2 h; (−8,98; −3,63; 1,125) lateral | 56,70 (224,2 h) | 1,54 | +1,5×10⁻⁵ | 4,1 min |
| 0,125×0,125×0,05 | 19,957 | **18,567** | 1027,8 h; (−8,98; −3,69; 0,025) base | 56,70 (224,4 h) | **1,43** | −2,7×10⁻⁵ | 17,6 min |
Mejora respecto de AI 50: −1,29 (0,25³), −1,38 (anisótropa), −1,39 °C (fina). Energía generada vs esperada +0,066 % en todas; volúmenes por vertida iguales a los del Anexo 5 (la aislación no los cambia). El sesgo de malla (fina − 0,25³) baja de +0,51 °C (AI 50) a +0,41 °C (AI 100). Con AI 100 el máximo se reparte entre una cara lateral y la base, por lo que aumentar más la aislación inferior rendiría poco.
**Precisión temporal.** Estar bajo el paso estable (1250 s anisótropa, 1023 s fina) garantiza estabilidad, no precisión. La comprobación temporal existente (dt 600 → 300 → 150 s; |ΔM_Δ| < 0,01 °C) se hizo con el plan de 63,7 m³/h, AI 50 y malla 0,25³: se aplica **por analogía**, no está medida para este candidato. Pendiente opcional (~1 min): AI 100, malla 0,25³, dt 300 y 150 s.
**Conclusión (un escenario crítico, suelo 0,8 m).** Con AI 100 mm el candidato queda por debajo de 20 °C en las tres mallas, con 1,43 °C de margen en la más fina, mayor que el sesgo de malla observado (0,41–0,51 °C). Mejora respecto de la versión AI 50 (margen 0,043 °C). **No reemplaza la evaluación completa de 24 escenarios × 2 suelos ni el retiro** en la configuración AI 100, no ejecutados. **Los 100 mm son una propuesta térmica condicionada: no se afirma que sea una solución ejecutable.** Falta justificar el material bajo la fundación (tipo de aislante, resistencia a compresión y fluencia frente a la presión de apoyo, drenaje, interacción con el hormigón de limpieza y la capacidad geotécnica) y su k.
**Costo actualizado de las campañas pendientes (con muestreo reducido, malla 0,25³, progresivo 36 m³/h, AI 100):** evaluación completa 24 escenarios × 2 suelos a 56 d: 48 corridas, ~0,5–1 min cada una → ~5–8 min con 12 procesos; retiro (una capa de 20 mm cada 14 d desde el día 56) hasta 150 d: 48 corridas, ~1,5 min cada una → ~10–15 min. No lanzadas.

---
# Anexo 7 — candidato con aislación inferior de 100 mm: campañas completas (18-11-2026)
Candidato: progresivo 36 m³/h (D = 17,69 h), 10 °C, PUR 80 mm (15 h), aislación inferior **100 mm**, mallas 0,25³, dt 600 s, muestreo reducido. 96 corridas (24 escenarios × suelos 2 m y 0,8 m × {56 d sin retiro, 150 d con retiro de capas de 20 mm en los días 56, 70, 84, 98}), ~30 min con 12 procesos, sin errores. Tablas: `ai100_planes_56d.csv`, `ai100_escenarios_56d.csv`, `ai100_retiro_planes_150d.csv`, `ai100_retiro_escenarios_150d.csv` (script `defensa_analisis_ai100.py`). **Supuesto climático:** solo hasta el día 14 hay pronóstico; 14–28 d trayectorias GEFS sin habilidad; **desde 28 d (hora 672) se repite el último ciclo diario**; todo el retiro (días 56–98) y los resultados a 150 d dependen de él.

**1. Precisión temporal (candidato actual, suelo 0,8 m, M4_SS cálida, 0,25³, 56 d).** dt = 600 / 300 / 150 s → ΔT máx 18,1572 / 18,1492 / 18,1491 °C. dt = 600 s sobreestima ΔT en 0,008 °C (conservador); 300→150: 0,00006 °C. Instante 1051,5 h en todos. Las comprobaciones previas (< 0,01 °C) son por tanto aplicables a este candidato. **Diferencia de 1,3×10⁻⁵ °C por muestreo reducido:** `simular_vertidas` integra por tramos entre eventos consecutivos (activaciones, cambios de protección *y los instantes de muestreo*) y subdivide cada tramo en pasos iguales ≤ dt; al quitar las salidas horarias los tramos cambian de longitud y con ellos el reparto de pasos, es decir los pasos de integración **sí** cambian ligeramente. La diferencia (1,3×10⁻⁵ °C) es del orden del error de truncamiento temporal y 600 veces menor que el efecto dt 600→300.

**2. 24 escenarios hasta 56 d (sin retiro), AI 100 mm:** P_T = P_Δ = P_f = 0 en las 3 interpretaciones, suelos 2 m y 0,8 m, a 14, 28 y 56 d.
| Suelo | Interpretación | máx ΔT 14 d / 28 d / 56 d [°C] | margen a 20 °C (56 d) | Tmax máx |
|---|---|---|---:|---:|
| 2 m | M1_SS | 13,08 / 14,23 / 15,43 | 4,57 | 51,06 |
| 2 m | 4 SS · 8 curvas | 15,19–16,41 / 16,71 / **18,01** (M4_SS cálida, 1123 h) | 1,99 | 56,61 |
| 0,8 m | M1_SS | 13,08 / 14,55 / 15,54 | 4,46 | 51,04 |
| 0,8 m | 4 SS · 8 curvas | 15,19–16,41 / 17,11 / **18,15** (M4_SS cálida, 1075,5 h) | 1,85 | 56,58 |
Frente al candidato AI 50 (Anexo 2): ΔT máx 18,42/19,45 → 18,01/18,15 °C (suelo 2 m / 0,8 m); en malla 0,25³ el escenario crítico sigue siendo M4_SS cálida.

**3. Retiro de capas (150 d), AI 100 mm: aquí cambia el resultado.**
| Suelo | Interpretación | máx ΔT 56 d → 150 d | P_Δ = P_f a 150 d | margen a 20 °C (150 d) |
|---|---|---|---|---:|
| 2 m | M1_SS | 15,43 → 18,15 | 0 | 1,85 |
| 2 m | 4 SS | 18,01 → **20,71** | **0,081** (> 5 %) | **−0,71** |
| 2 m | 8 curvas | 18,01 → **20,71** | **0,040** (un escenario falla: M4_SS cálida, peso 0,0403; formalmente ≤ 5 %) | **−0,71** |
| 0,8 m | M1_SS · 4 SS · 8 curvas | 15,54/18,15 → 17,16 / 19,56 / 19,56 | 0 | 2,84 / 0,44 / 0,44 |
- **El escenario crítico cambia:** sin retiro es M4_SS cálida, en la base (día ~43–47); con retiro es M4_SS cálida en la **cara lateral x≈+8,98, y≈+3,63, z≈1,125 m, a las 2394 h (día 99,8), 2,4 días después de retirar la cuarta capa**, en los 24 escenarios (todos superan su máximo de 56 d; exceso máx 2,72 °C con suelo 2 m y 1,62 °C con 0,8 m). Máximos tras cada retiro, peor escenario (suelo 2 m · 0,8 m): 18,29 · 18,17; 18,55 · 18,13; 18,88 · 18,11; **20,69 · 19,54 °C**.
- Con AI 50 el mismo retiro daba 18,73 °C (suelo 2 m) y 16,77 °C tras el retiro completo (suelo 0,8 m): **más aislación inferior empeora el pico posterior al retiro** (el núcleo, menos enfriado por la base, queda más caliente cuando se expone la cara lateral). Ni AI 50 ni AI 100 cumplen holgadamente con ambos suelos bajo esta regla: AI 50 falla con 0,8 m en malla fina (0,043 °C de margen, sin retiro); AI 100 falla con suelo 2 m con retiro (−0,71 °C) y deja 0,44 °C con 0,8 m.
- Tendencia final (días 140–150): ΔT decrece (−0,07 a −0,17 °C/día); ΔT final ≤ 8,4 °C; Tmax del núcleo 21–26 °C. No hay crecimiento al final; el máximo es el pico posterior al último retiro.
- M_T máx 56,6 °C en todos los casos (D = 17,69 h).

**4. Escenario más exigente que el verificado en malla fina.** Sí: el pico lateral tras el último retiro (20,71 / 19,56 °C) **no** se comprobó en malla fina (la malla fina se verificó para el pico en la base sin retiro). Puede cambiar la decisión (suelo 2 m ya supera 20 °C en malla 0,25³, y para 0,8 m el margen 0,44 °C es menor que el sesgo de malla 0,4–0,5 °C). Comprobación dirigida propuesta (no lanzada): (a) M4_SS cálida, AI 100, suelo 0,8 m, con retiro hasta ~3000 h, en 0,25×0,25×0,05 (~11 min) y en 0,125×0,125×0,05 (~45 min, ~4 GB); (b) ajuste de la regla de retiro, medido en malla 0,25³ y los dos suelos en M4_SS cálida: capas retiradas cada 21 d o 28 d, o conservar la última capa hasta el día ≥ 120 (4–6 variantes × ~1,5 min ≈ 6–10 min); (c) si alguna variante deja ≥ 1,5 °C de margen en ambos suelos, repetir la matriz de 24 escenarios para esa regla (~10–15 min). Total estimado ≈ 30–40 min + ~45 min de malla fina.

**5. Qué está cerrado y qué falta (para el jueves)**
- *Cumplimiento en los escenarios calculados (malla 0,25³):* candidato 36 m³/h progresivo con 10 °C y PUR 80 mm: **sin retiro** cumple (P_f = 0) con AI 50 o AI 100 a 14/28/56 d en ambos suelos, bajo las 3 interpretaciones. **Con la regla de retiro probada** (una capa de 20 mm cada 14 d desde el día 56): cumple con AI 50 (suelo 2 m, 150 d) pero **no** con AI 100 en suelo 2 m (P_f = 0,04–0,08, ΔT 20,71 °C) y con margen 0,44 °C en suelo 0,8 m.
- *Sensibilidad numérica:* temporal ≈ 0,008 °C (dt 600 conservador); espacial +0,4–0,5 °C en el pico de la base (AI 50: 19,45 → 19,96 °C, suelo 0,8 m; AI 100: 18,16 → 18,57 °C); el pico lateral post-retiro no tiene verificación espacial.
- *Pendiente constructivo:* viabilidad estructural y material de la aislación inferior (50/100 mm), enfriamiento a 10 °C, caudal de 36 m³/h (el modelo activa 31–42 m³/h numéricos con el mismo cronograma), producto y k de la manta, regla de retiro ejecutable, fundamento del suelo (2 m / 0,8 m), y significado de las 8 curvas.
- *Alcance de P_f = 0:* conjunto discreto de 24 escenarios (3 trayectorias GEFS × 8 curvas), pesos declarados, malla 0,25 m, extensión climática supuesta tras el pronóstico. No es un resultado de fisuración.

---
# Anexo 8 — reglas de retiro más gradual o tardío (aislación inferior 100 mm; M4_SS cálida; ambos suelos)
Candidato idéntico al Anexo 7 (progresivo 36 m³/h, 10 °C, PUR 80 mm en 4 capas de 20 mm, AI 100 mm, malla 0,25³, dt 600 s, muestreo reducido). 14 corridas (7 reglas × 2 suelos), horizonte 5400 h (225 d) para capturar el pico tras el retiro completo, ~11 min con 14 procesos, sin errores. Tabla `ret2_reglas_M4SS_calida.csv`. **Supuesto climático:** desde la hora 672 se repite el último ciclo diario; todo lo posterior al día 28 depende de él. Fechas desde t0 = jueves 19-11-2026 10:00 (Ovalle).
| Regla | Días de retiro (capa 1–4) | Fechas | Protección total | Máx ΔT suelo 2 m (día, ubicación) | Máx ΔT suelo 0,8 m | Margen 2 m / 0,8 m |
|---|---|---|---:|---|---|---|
| **A (anterior, falla)**: 14 d desde d56 | 56, 70, 84, 98 | 14-01, 28-01, 11-02, 25-02-2027 | 98 d | **20,71 (d99,8, lateral)** | 19,56 (d99,8, lateral) | **−0,71** / 0,44 |
| B: 21 d desde d56 | 56, 77, 98, 119 | 14-01, 04-02, 25-02, 18-03 | 119 d | 18,29 (d60,8, lateral) | 18,17 (d59,8) | 1,71 / 1,83 |
| C: 28 d desde d56 | 56, 84, 112, 140 | 14-01, 11-02, 11-03, 08-04 | 140 d | 18,29 (d60,8) | 18,17 (d59,8) | 1,71 / 1,83 |
| D: 14 d, última en d140 | 56, 70, 84, 140 | 14-01, 28-01, 11-02, 08-04 | 140 d | 18,89 (d87,8) | 18,17 (d59,8) | 1,11 / 1,83 |
| E: 14 d, última en d168 | 56, 70, 84, 168 | …, 06-05 | 168 d | 18,89 (d87,8) | 18,17 (d59,8) | 1,11 / 1,83 |
| F: 21 d desde d84 | 84, 105, 126, 147 | 11-02, 04-03, 25-03, 15-04 | 147 d | 18,01 (d45,8, base sin retiro) | 18,16 (d43,8) | **1,99 / 1,84** |
| G: 28 d, última en d168 | 56, 84, 112, 168 | 14-01, 11-02, 11-03, 06-05 | 168 d | 18,29 (d60,8) | 18,17 | 1,71 / 1,83 |
Ubicación del pico lateral: cara x ≈ +8,98 (o −8,98), y ≈ +3,63 (o −3,63), z = 1,125 m. Máximos tras el retiro completo (suelo 2 m · 0,8 m): A 20,69 · 19,54; B 16,61 · 15,48; C 13,25 · 12,22; D 12,02 · 11,06; E 8,31 · 7,63; F 13,31 · 12,27; G 9,23 · 8,46. En todos ΔT decrece al final (−0,01 a −0,04 °C/d; ΔT final ≤ 4,0 °C a 225 d): el horizonte captura el pico.
**Lectura.** (1) La regla A falla (ΔT 20,71 °C, suelo 2 m) porque la **cuarta capa se retira solo 14 d después de la tercera**; separar las retiradas 21 d (B) o 28 d (C) baja el pico posterior al retiro completo a 16,6/13,3 °C (2 m). (2) En B–E y G el máximo pasa a ser el pico del **primer retiro (día 56)**, común a todas (18,29 / 18,17 °C), y la regla D/E (retrasar la última capa) mantiene el pico de la tercera retirada en 18,89 °C con suelo 2 m (margen 1,11 °C, **< 1,5 °C**). (3) La regla F (retiro desde el día 84) no genera un pico nuevo: el máximo es el pico en la base sin retiro (18,0–18,2 °C) que ya se verificó en malla fina con AI 100 (0,125×0,125×0,05: 18,57 °C, margen 1,43 °C).
**Reglas prometedoras (margen ≥ 1,5 °C en ambos suelos, objetivo de selección y no garantía):** **B** (protección más corta: 119 d, margen 1,71/1,83), **C** (140 d, 1,71/1,83) y **F** (147 d, 1,99/1,84; retira más tarde, sin pico nuevo). D, E (1,11 °C en suelo 2 m) y A quedan descartadas por esta selección; G equivale a C pero con 168 d de protección.
**Costo de verificar el nuevo pico crítico (no lanzado):** para B y C el pico es el del primer retiro (cara lateral, días ~60), idéntico en ambas, y las capas posteriores lo reducen: basta una corrida por suelo hasta ~d125 (3000 h; cubre el pico a d60 y el retiro final de B) en cada malla. Malla 0,25×0,25×0,05: ~9 min por corrida (suelo 0,8 m; la de 2 m, ~9 min); malla 0,125×0,125×0,05: ~40 min y ~4 GB por corrida. Propuesta: B con suelo 0,8 m en ambas mallas (≈ 50 min) y, si se prefiere el margen, F con suelo 0,8 m en la anisótropa para el pico post-retiro (~11 min, al ser este mucho menor que el pico base ya verificado). Después, evaluar los 24 escenarios × 2 suelos con la regla escogida (~25–30 min a 3400 h con 0,25³). Ninguna de estas campañas se ha lanzado.

---
# Anexo 9 — verificación espacial de la regla B (calendario de retiro cada 21 d desde el día 56)
**Naturaleza de la regla.** B es un **calendario evaluado** (capas de 20 mm retiradas el 14-01, 04-02, 25-02 y 18-03-2027; días 56, 77, 98, 119 desde t0 = 19-11-2026 10:00), simulado bajo la **extensión climática supuesta** (pronóstico hasta el día 14; trayectorias GEFS sin habilidad hasta el 28; ciclo diario repetido desde la hora 672). **La retirada real requerirá seguimiento térmico y un criterio de autorización; no se ha validado ninguna regla basada en mediciones.** F (retiros días 84, 105, 126, 147) se mantiene como alternativa **sin verificar** (sus retiros no se han comprobado en mallas finas; solo se verificó el pico base sin retiro).
Candidato: progresivo 36 m³/h, 10 °C, PUR 80 mm (4 capas), AI 100 mm, M4_SS cálida, horizonte 3400 h (141,7 d), dt 600 s, muestreo reducido, métricas en cada paso. Corridas: anisótropa 0,25×0,25×0,05 (ambos suelos, 10,8 min c/u en paralelo) y fina 0,125×0,125×0,05 (suelo 2 m, 38,9 min, ~4 GB). Tabla comparativa con la malla gruesa (Anexo 8).
| Suelo / malla | Máx ΔT global (día; ubicación S25) | Máx tras retiro 1 / 2 / 3 / 4 (completo) | Margen a 20 °C | Tmax |
|---|---|---|---:|---:|
| 2 m · 0,25³ | 18,29 (d60,8; lateral x=+8,98 y=+3,63 z=1,125) | 18,29 / 17,78 / 16,74 / 16,61 | 1,71 | 56,6 |
| 2 m · aniso | 18,58 (d60,8; −8,98, −3,63, 1,125) | 18,58 / 18,09 / 17,06 / 17,26 | 1,42 | 56,72 |
| **2 m · fina** | **18,63 (d60,8; −8,98, −3,69, 1,075)** | 18,63 / 18,12 / 17,09 / 17,38 | **1,37** | 56,72 (d9,4) |
| 0,8 m · 0,25³ | 18,17 (d59,8; lateral) | 18,17 / 17,29 / 15,88 / 15,48 | 1,83 | 56,6 |
| 0,8 m · aniso | 18,47 (d59,8; −8,98, −3,63, 1,125) | 18,47 / 17,60 / 16,20 / 16,12 | 1,53 | 56,70 |
Pre-retiro (días 0–56): 18,37 (2 m, fina), 18,46 (0,8 m, aniso). Balance de energía: residuo ≤ 3,8×10⁻⁵ J en todas; energía generada +0,066 % vs esperada. Paso estable 1250 s (aniso) y 1023 s (fina) > 600 s; la sensibilidad temporal medida es 0,008 °C (dt 600 conservador). ΔT al final (d141,7) = 11–12 °C y todavía decreciente (−0,22 °C/d): el pico quedó capturado.
**Lectura.** (1) Con suelo 2 m el pico del primer retiro crece +0,29 °C (aniso) y +0,34 °C (fina) sobre la malla gruesa, y se estabiliza entre las dos mallas finas (+0,05 °C); margen 1,37 °C, mayor que la sensibilidad numérica observada (0,3–0,5 °C), pero no es una cota. (2) Los retiros 2–4 quedan 0,5–1,5 °C por debajo del primero en todas las mallas: **no pasan a controlar**. (3) Con suelo 0,8 m el pico del primer retiro (18,47 °C aniso) es comparable al pico base sin retiro verificado en malla fina (18,57 °C): ambos con margen ≈ 1,4–1,5 °C; con el incremento aniso→fina observado en suelo 2 m (+0,05 °C) y en el pico base (+0,11 °C) se estima 18,5–18,6 °C en malla fina (estimación, no corrida). No hay ningún pico cercano a 20 °C que justifique otra corrida fina. **B conserva margen frente a la sensibilidad numérica observada.**
**Siguiente paso propuesto (no lanzado):** matriz completa de 24 escenarios × suelos 2 m y 0,8 m con la regla B, AI 100 mm, malla 0,25³ y horizonte 3400 h; 48 corridas con muestreo reducido, ~15–18 min con 12 procesos. Sus resultados deben leerse con un sesgo de malla de +0,3–0,5 °C en el ΔT; no sustituyen la verificación fina del escenario crítico (hecha solo con suelo 2 m para el pico de retiro y 0,8 m para el pico base).

---
# Anexo 10 — regla B: matriz completa de 24 escenarios × 2 suelos (18-11-2026)
Candidato exacto: colocación continua **progresiva a 36 m³/h** (D = 17,69 h), 10 °C, manta PUR de 80 mm en 4 capas de 20 mm (desde las 15 h), aislación inferior 100 mm (k = 0,04 supuesto), **retiro de 20 mm los días 56, 77, 98 y 119** (14-01, 04-02, 25-02, 18-03-2027), malla 0,25³, dt 600 s, muestreo reducido, horizonte 3400 h (141,7 d). 48 corridas, ~19 min con 12 procesos, sin errores. Tablas: `retB_escenarios.csv` (por escenario: D, M_T, ΔT máx, instantes, ubicaciones y margen a 14 d, 28 d, 56 d y total; máximos tras cada retiro; tendencia final) y `retB_planes.csv` (P_T, P_Δ y P_f por interpretación y horizonte). Script `defensa_analisis_retB.py`. **Calendario evaluado bajo la extensión climática supuesta** (pronóstico hasta el día 14; trayectorias GEFS sin habilidad hasta el 28; ciclo diario repetido desde la hora 672); la retirada real requiere seguimiento térmico y un criterio de autorización aún no validados.
## Resultados calculados (malla 0,25³)
| Suelo | Interpretación | P_T | P_Δ | P_f (14 d / 28 d / 56 d / 141,7 d) | Máx ΔT 56 d → 141,7 d [°C] | Margen a 20 °C | Máx M_T [°C] |
|---|---|---|---|---|---|---:|---:|
| 2 m | M1_SS | 0 | 0 | 0 / 0 / 0 / 0 | 15,43 → 15,71 | 4,29 | 51,06 |
| 2 m | 4 SS | 0 | 0 | 0 / 0 / 0 / 0 | 18,01 → 18,29 | 1,71 | 56,61 |
| 2 m | 8 curvas | 0 | 0 | 0 / 0 / 0 / 0 | 18,01 → 18,29 | 1,71 | 56,61 |
| 0,8 m | M1_SS | 0 | 0 | 0 / 0 / 0 / 0 | 15,54 → 15,60 | 4,40 | 51,04 |
| 0,8 m | 4 SS · 8 curvas | 0 | 0 | 0 / 0 / 0 / 0 | 18,15 → 18,17 | 1,83 | 56,58 |
P_f = 0 corresponde al conjunto discreto evaluado (3 trayectorias GEFS × 8 curvas, pesos declarados), no a una probabilidad real de fisuración. Ningún escenario falla (0/24 por suelo).
**Escenario crítico.** M4_SS cálida en los dos suelos; es el mismo escenario verificado en malla fina. Suelo 2 m: máx 18,29 °C el día 60,8, cara lateral (x = +8,97, y = +3,62, z = 1,12 m), tras el primer retiro (controla en 12 de 24 escenarios; en los otros 12 controla el pico en la base sin retiro, 18,01 °C). Suelo 0,8 m: máx 18,17 °C el día 59,8 (primer retiro, lateral), seguido del pico en la base (18,15 °C, día ~43; controla en 20 de 24 escenarios, p. ej. M4_SS central y fría a los 40,8 d). Picos tras cada retiro (máximo sobre los 24 escenarios, 2 m · 0,8 m): 18,29 · 18,17; 17,78 · 17,29; 16,74 · 15,88; 16,61 · 15,48. **El escenario y la ubicación crítica no cambian respecto de los verificados en malla fina (Anexo 9), y no aparece ningún caso más exigente.** M_T máx 56,6 °C (día ~9,3, interior); D no cambia (17,69 h).
**Tendencia final y horizonte.** Ningún pico cae en los últimos 10 d del horizonte; instantes de máximo entre el día 1,1 y el 60,8. A 141,7 d el ΔT decrece en todos los escenarios (pendiente −0,11 a −0,22 °C/d; ΔT final ≤ 11,9 °C), de modo que el pico quedó dentro del horizonte; si los retiros se retrasaran, habría que ampliar el horizonte.
## Sensibilidad espacial observada (separada de lo calculado)
Los valores de arriba son de la malla 0,25³. En el escenario crítico se **observó** al refinar (B, M4_SS cálida): suelo 2 m: 18,29 (0,25³) → 18,58 (anisótropa) → 18,63 °C (0,125×0,125×0,05; margen 1,37 °C); suelo 0,8 m: 18,17 → 18,47 °C (anisótropa; margen 1,53 °C), con pico base fino ya verificado en 18,57 °C (Anexo 6). Estas diferencias (+0,3 a +0,5 °C) son observaciones en un escenario, **no se han aplicado ni extrapolado a los otros 23 escenarios y no son una cota ni una corrección conocida**. Los márgenes calculados en malla gruesa (1,71/1,83 °C) y los observados en el escenario crítico más fino (1,37/1,53 °C) se informan por separado.
## Qué está cerrado y qué falta para el jueves
**Cerrado (en los escenarios calculados):** candidato B (36 m³/h progresivo, 10 °C, PUR 80 mm, AI 100 mm, retiro días 56/77/98/119): sin fallos en 24 escenarios × 2 suelos, P_f = 0 en las tres interpretaciones de pesos; escenario crítico verificado en malla anisótropa (ambos suelos) y fina (suelo 2 m para el pico de retiro; suelo 0,8 m para el pico base); sensibilidad temporal 0,008 °C; balance de energía cerrado (≤ 4×10⁻⁵ J); la regla anterior (retiros cada 14 d: ΔT 20,71 °C) está documentada como fallida; el desempate del voxelizador no se aplicó.
**Pendiente (no resuelto por simulación):** (1) viabilidad estructural y justificación del material de 100 mm bajo la fundación (compresión, fluencia, drenaje, k); (2) enfriamiento a 10 °C; (3) caudal real de 36 m³/h (el modelo activa 31–42 m³/h numéricos); (4) producto y k de la manta, y su retiro por capas en obra; (5) criterio de autorización de retiro basado en mediciones (el calendario es solo evaluado); (6) fundamento del suelo efectivo (2 m / 0,8 m, T = 20 °C); (7) significado de las 8 curvas y la mezcla del IFC; (8) pronóstico solo hasta el día 14: lo posterior es supuesto; (9) F sin verificar; (10) no hay verificación fina del pico de retiro con suelo 0,8 m ni de los otros escenarios.
## Síntesis reproducible (PowerShell, desde la raíz, `$env:PYTHONPATH="src"; $env:PYTHONIOENCODING="utf-8"`)
```
python scripts/preparar_clima_extendido_gefs.py
python scripts/defensa_comprobaciones.py --procs 12 retB            # 48 corridas, reanudable
python scripts/defensa_analisis_retB.py                             # tablas retB_*.csv
python scripts/defensa_malla_candidato.py correr aniso005 100 reducido 2.0 B 3400
python scripts/defensa_malla_candidato.py correr comb 100 reducido 2.0 B 3400   # ~39 min, ~4 GB
python scripts/defensa_malla_candidato.py correr comb 100 reducido 0.8 none 1344 # pico base fino, suelo 0,8 m
python -m unittest discover -s tests                                 # 189 pruebas
```
Sin commit ni push. Archivos nuevos o modificados sin versionar: `scripts/defensa_*.py`, `scripts/preparar_clima_extendido_gefs.py`, `scripts/auditar_curvas_hidratacion.py`, `tests/test_defensa.py`, `configs/defensa/`, `reports/Defensa/`, `data/weather/…_28d/`, y cambios en `src/mcplan/thermal/decision.py`.
