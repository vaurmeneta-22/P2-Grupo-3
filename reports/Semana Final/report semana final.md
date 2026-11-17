# Report semana final — trabajo y simulaciones del 16 y 17 de noviembre de 2026

## Alcance

Este informe resume las decisiones, verificaciones y corridas térmicas realizadas desde el 16 de noviembre para evaluar el plan de hormigonado de la fundación. Las simulaciones finales consideradas aquí usan el IFC del proyecto, un horizonte de 14 días y Δx = 0,25 m. Los resultados de exploración que preceden a la corrección del seguimiento S25 se identifican como preliminares.

Presentación navegable: [abrir TermoBIM · Fundación Grupo 3](presentacion_final.html).

## Configuración constructiva evaluada

- Inicio de obra: jueves 19 de noviembre de 2026, 10:00, hora de Ovalle (UTC−3).
- Secuencia de 10 h: V1 0–2 h; V2 2–4 h; V3 4–7 h; V4 7–9 h; V5 9–10 h. Es un caudal medio de aproximadamente 63,7 m³/h para los 636,8 m³ del IFC.
- Temperatura de colocación en las corridas recientes: 10 °C en todas las vertidas. En una exploración anterior se había propuesto 15 °C como referencia tras comparar temperaturas de colocación; 10 °C fue una variante posterior para buscar reducir el diferencial.
- Hormigón: propiedades térmicas M1, k = 2,0 W/(m·K), ρ = 2400 kg/m³ y cp = 900 J/(kg·K). La fuente nominal es M1_SS, tomada de `data/hidratacion/source_terms.txt`.
- Malla: voxeles de 0,25 m; paso temporal máximo 600 s; duración de 336 h.
- Clima: snapshots GEFS seleccionados como frío, central y cálido. No se simularon los 31 miembros individuales del ensamble.
- Intercambio con el aire: Robin, h = 8 W/(m²·K).
- Base: suelo efectivo de 2 m, k = 1,5 W/(m·K) y temperatura profunda constante de 20 °C; se representa como condición equivalente y no como voxeles de suelo.
- Moldaje lateral: V1, V2 y V5, 18 mm, k = 0,1 W/(m·K), retiro según las 72 h del IFC desde la colocación. V3 y V4 sin moldaje.
- Manta PUR: k = 0,025 W/(m·K), sobre caras laterales y superiores, excluyendo la cara inferior. Se instala a t = 15 h. En la referencia de 40 mm no se retira dentro del horizonte de 14 días.

## Verificaciones y ajustes técnicos

Se detectó que el observador S25 podía conservar superficies de una vertida o protección que ya no correspondían al estado actual. Se corrigió la clave de caché para considerar la máscara activa y las resistencias de borde, y se recalcularon los casos afectados. La verificación registrada informa 184 pruebas automatizadas aprobadas.

También se revisó la ejecución de E11 en Windows: la consola podía fallar al imprimir símbolos por codificación. La salida UTF-8 permitió ejecutar y guardar el resumen; no era un error del solver. Esa auditoría reportó Tmax = 38,679 °C en su propio caso E11, que usa otro horizonte y no forma parte de la matriz de la fundación.

Se probaron variantes de hora de inicio (10:00 y 22:00), temperatura de colocación, profundidad efectiva del suelo y espesor de manta. Las primeras tablas de ΔT que usaban la caché S25 anterior son exploratorias y no deben usarse para decidir cumplimiento; los resultados centrales de espesores se recalcularon con el seguimiento corregido. La malla final de estas comparaciones es 0,25 m; Δx = 0,063 m quedó descartado.

## Comparación central M1_SS / clima GEFS central

Las cuatro corridas mantienen el mismo plan, clima y demás parámetros; solo cambia el espesor PUR. El ΔT es el máximo global Tmax del hormigón menos la temperatura S25 mínima reconstruida en la superficie, cuyos puntos pueden estar separados espacialmente.

| Manta PUR | Tmax [°C] | Hora de Tmax [h] | ΔT máximo [°C] | Hora del ΔT máximo [h] | Evaluación individual |
|---:|---:|---:|---:|---:|---|
| 30 mm | 50,692 | 207,17 | 19,749 | 333,5 | Cumple por 0,251 °C |
| 40 mm | 50,735 | 208,83 | 17,912 | 336 | Cumple dentro del horizonte |
| 50 mm | 50,774 | 210,33 | 16,730 | 336 | Cumple dentro del horizonte |
| 60 mm | 50,802 | 211,67 | 15,832 | 336 | Cumple dentro del horizonte |

La tendencia del caso central es clara: más aislación reduce ΔT y sube levemente Tmax. Todas estas Tmax quedan muy por debajo de 70 °C. En las mantas de 40, 50 y 60 mm, el máximo de ΔT aparece al final del horizonte, por lo que no está demostrado que siga bajo el límite después del día 14.

### Controles adicionales de manta de 40 mm

| Variante | ΔT máximo [°C] | Instante [h] | Tmax [°C] | Resultado |
|---|---:|---:|---:|---|
| Sin manta PUR | 32,567 | 164,0 | 50,363 | Excede ΔT de 20 °C |
| PUR 40 mm, retirada a las 240 h | 31,784 | 330,17 | 50,735 | Excede ΔT de 20 °C; variante archivada, no adoptada |
| PUR 40 mm, sin retiro durante los 14 días | 17,912 | 336,0 | 50,735 | Caso de referencia central |

La variante de retiro al día 10 se conserva como experimento, pero se descartó como configuración viable. La corrida sin manta también muestra que la temperatura máxima no es el problema principal: lo es la diferencia térmica global.

## Matriz completa: manta PUR de 40 mm

Se corrieron 24 combinaciones para el plan de 10 h, colocación a 10 °C y manta de 40 mm. Se asignó peso uniforme a las ocho columnas de hidratación y se usaron los pesos configurados para GEFS frío, central y cálido. El límite evaluado es ΔT ≤ 20 °C y Tmax ≤ 70 °C; el criterio agregado de fallo es P_f ≤ 5 %.

- 20 de 24 combinaciones cumplen ambos límites de forma individual.
- Las 24 cumplen Tmax ≤ 70 °C; la mayor Tmax de la matriz es aproximadamente 56,21 °C.
- Cuatro exceden ΔT: M2_SS con clima frío (20,31 °C) y M4_SS con clima frío, central y cálido (21,45; 21,13; 20,80 °C).
- Ponderando los escenarios, P_T = 0 %, P_ΔT = P_f = 16,94 %. Por tanto, este plan no cumple el criterio agregado P_f ≤ 5 %.

La matriz incluye M1–M4 de las campañas SS y 22A. Como el IFC identifica la fundación como M1, M1_SS es el caso nominal más directamente alineado; las curvas M2–M4 deben interpretarse como sensibilidad de hidratación, no como mezclas confirmadas de esta fundación. M1_22A también es una curva de sensibilidad y tiene menor confianza indicada en el archivo fuente.

### Repetición completa con manta PUR de 60 mm

Se repitieron las mismas 24 combinaciones y el mismo plan, cambiando únicamente el espesor de la manta de 40 a 60 mm.

- **24 de 24 casos cumplen** ΔT ≤ 20 °C y Tmax ≤ 70 °C.
- P_T = 0 %, P_ΔT = 0 % y P_f = 0 %, menor que el máximo permitido de 5 %.
- El peor ΔT es **18,964 °C**, en M4_SS con clima frío, a las 336 h.
- La mayor Tmax es **56,298 °C**, en M4_SS con clima central, a las 214 h.
- ΔT disminuyó en los 24 casos respecto de 40 mm, con una reducción media de **1,95 °C**. Los cuatro casos que fallaban con 40 mm pasan a cumplir con 60 mm.

Así, 60 mm es el primer espesor evaluado que cumple el criterio agregado para las ocho curvas y los tres climas. El máximo diferencial del peor caso ocurre en la última hora simulada; el resultado no demuestra el comportamiento posterior a los 14 días. La matriz también incluye curvas M2–M4 como sensibilidad, aunque el IFC clasifica la fundación como M1.

## Comparación de planes de colocación

Se evaluaron seis alternativas con la misma geometría IFC, propiedades térmicas, suelo, límite de 100 m³/h por cuadrilla y ventana climática. Cada plan se cruzó con las ocho curvas de hidratación y las tres trayectorias GEFS, para un total de **144 simulaciones (24 por plan)**. El inicio común es el 19 de noviembre de 2026 a las 10:00, hora local de Ovalle. Se aplican los límites ΔT ≤ 20 °C, Tmax ≤ 70 °C y Pf ≤ 5 %.

| Plan | Duración D [h] | Tmax peor caso [°C] | ΔT peor caso [°C] | Pf | Resultado |
| --- | ---: | ---: | ---: | ---: | --- |
| Secuencia original IFC, 18 °C, sin manta | 24,384 | 63,269 | 45,886 | 100 % | No cumple |
| Programa de 10 h, 10 °C, sin manta | 10,000 | 55,740 | 38,662 | 100 % | No cumple |
| Programa de 10 h, 10 °C, PUR 60 mm | 10,000 | 56,298 | 18,964 | 0 % | Cumple |
| Secuencia continua a 100 m³/h, PUR 60 mm | 6,368 | 56,253 | 18,899 | 0 % | Cumple |
| Secuencia continua de 12 h, PUR 60 mm | 12,000 | 56,363 | 19,003 | 0 % | Cumple |
| Activación conjunta idealizada, PUR 60 mm | 6,368 | 56,115 | 18,834 | 0 % | Cumple en el modelo |

“Peor caso” toma el máximo de cada métrica entre los 24 escenarios; por ello Tmax y ΔT pueden corresponder a curvas o climas distintos. Los cuatro planes con manta cumplen en las 24 combinaciones: Tmax queda bajo 70 °C y ΔT bajo 20 °C; sus máximos de ΔT ocurren en el escenario frío M4_SS al final de la simulación (336 h). Los dos controles sin manta fallan por ΔT en todos los escenarios. Para el plan de 10 h, la comparación con/sin manta cambia solo la protección y muestra el efecto de PUR de 60 mm bajo el mismo cronograma y temperatura de colocación.

Con 100 m³/h como capacidad máxima, el límite aritmético inferior para colocar 636,775 m³ es **6,368 h**. La secuencia continua sin esperas alcanza ese valor en los planes probados y cumple el criterio térmico; se considera la alternativa más rápida factible **entre las evaluadas**, no una prueba de óptimo global. La activación conjunta de las cinco regiones IFC también da D=6,368 h para comparar el límite inferior, pero en el cálculo térmico activa todo el volumen en t=0. Es una idealización y no una secuencia operativa realizable; no se la presenta como plan recomendado. Las duraciones describen el programa y el cálculo de D: cada vertida se activa instantáneamente al inicio de su tramo, simplificación ya adoptada para este solver.

La secuencia continua programa cada vertida sin espera, con duración proporcional a su volumen: V1=1,476 h; V2=1,476 h; V3=2,142 h; V4=0,889 h y V5=0,384 h. El plan de 12 h conserva la secuencia, a menor caudal medio (~53,1 m³/h). En todos los planes con manta, esta se instala cinco horas después de terminar la colocación programada; se mantienen el moldaje lateral de V1, V2 y V5, el suelo efectivo de 2 m a 20 °C y PUR de 60 mm (k=0,025 W/(m·K)) en las caras expuestas salvo la base.

La repetición de la matriz de referencia del plan 10 h / PUR 60 mm reprodujo exactamente sus métricas anteriores (diferencia máxima numérica 0,0 °C en Tmax y ΔT). Las series temporales ilustrativas adicionales usan M1_SS y clima central; no sustituyen la matriz de 144 casos.

### Figuras de comparación

- [Resumen de duración, Pf, ΔT y Tmax](evidencia_semana_final/comparacion_100m3h/resumen_planes.png)
- [Matriz de ΔT para los 24 escenarios](evidencia_semana_final/comparacion_100m3h/matriz_diferencial_escenarios.png)
- [Cronogramas de colocación e instalación de manta](evidencia_semana_final/comparacion_100m3h/cronogramas_planes.png)
- [Evolución de temperatura y ΔT para M1_SS / clima central](evidencia_semana_final/comparacion_100m3h/evolucion_temperaturas_y_delta.png)
- [Cortes térmicos del plan de 60 mm](evidencia_semana_final/comparacion_100m3h/cortes_temperatura_60mm.png)
- [Curvas de hidratación y climas de entrada](evidencia_semana_final/comparacion_100m3h/curvas_hidratacion_y_clima.png)
- [Geometría IFC](evidencia_semana_final/comparacion_100m3h/geometria_fundacion_ifc.png) · [geometría voxelizada](evidencia_semana_final/comparacion_100m3h/geometria_voxelizada_por_vertida.png)
- [Series temporales representativas de los seis planes (CSV)](evidencia_semana_final/comparacion_100m3h/series_temporales_representativas.csv)

## Archivos de resultados

- [Comparación de los seis planes (CSV)](evidencia_semana_final/comparacion_100m3h/resultados/comparacion_planes.csv)
- [Probabilidades ponderadas por plan (CSV)](evidencia_semana_final/comparacion_100m3h/resultados/probabilidades_por_plan.csv)
- [Resultados de los 144 escenarios (CSV)](evidencia_semana_final/comparacion_100m3h/resultados/resultados_por_escenario.csv)
- [Manifiesto de la corrida y trazabilidad](evidencia_semana_final/comparacion_100m3h/resultados/manifest.json)
- [Configuración reproducible de planes y corrida](../../configs/defensa_final_planes_100m3h.json)

- [Matriz de 24 resultados por escenario](evidencia_semana_final/matriz_8_curvas_3_climas/resultados/resultados_por_escenario.csv)
- [Probabilidades ponderadas y cumplimiento agregado](evidencia_semana_final/matriz_8_curvas_3_climas/resultados/probabilidades_por_plan.csv)
- [Dossier con la configuración de la matriz](evidencia_semana_final/matriz_8_curvas_3_climas/dossier_corrida.json)
- [Serie térmica M1_SS / GEFS central con PUR 40 mm](evidencia_semana_final/manta_40mm/series_temperatura.csv)
- [Gráfico térmico del caso de referencia 40 mm](evidencia_semana_final/manta_40mm/evolucion_temperatura_40mm.png)
- [Ubicación de extremos del caso de referencia 40 mm](evidencia_semana_final/manta_40mm/ubicacion_extremos_40mm.png)
- [Comparación documentada de mantas de 30, 50 y 60 mm](evidencia_semana_final/documentacion/comparacion_mantas_30_50_60.md)
- [Informe del recálculo S25](evidencia_semana_final/documentacion/resumen_recalculo_s25.md)
- [Auditoría de salida E11 en UTF-8](evidencia_semana_final/auditoria_e11_utf8/resumen_e11.txt)
- [Matriz completa con manta de 60 mm](evidencia_semana_final/matriz_60mm_8curvas_3climas/resultados/resultados_por_escenario.csv)
- [Probabilidades ponderadas de la matriz de 60 mm](evidencia_semana_final/matriz_60mm_8curvas_3climas/resultados/probabilidades_por_plan.csv)
- [Dossier y manifiesto de la matriz de 60 mm](evidencia_semana_final/matriz_60mm_8curvas_3climas/dossier_corrida_60mm.json) · [manifest](evidencia_semana_final/matriz_60mm_8curvas_3climas/resultados/manifest.json)
- [Comparación 40 mm vs. 60 mm (CSV)](evidencia_semana_final/matriz_60mm_8curvas_3climas/comparacion_40mm_vs_60mm.csv)
- [Gráfico comparativo 40 mm vs. 60 mm](evidencia_semana_final/matriz_60mm_8curvas_3climas/comparacion_40mm_vs_60mm.png)

### Nota de ejecución

Las dos matrices completaron sus 24 simulaciones y se guardaron en CSV. La corrida de 40 mm no generó manifiesto porque el dossier estaba fuera de la raíz del repositorio; el fallo ocurrió al final, luego de exportar las métricas. La corrida de 60 mm sí tiene manifiesto y se incluye en la evidencia.
