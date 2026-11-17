# Report semana final — trabajo y simulaciones del 16 y 17 de noviembre de 2026

## Alcance

Este informe resume las decisiones, verificaciones y corridas térmicas realizadas desde el 16 de noviembre para evaluar el plan de hormigonado de la fundación. Las simulaciones finales consideradas aquí usan el IFC del proyecto, un horizonte de 14 días y Δx = 0,25 m. Los resultados de exploración que preceden a la corrección del seguimiento S25 se identifican como preliminares.

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

## Matriz completa: ocho curvas × tres climas

Se corrieron 24 combinaciones para el plan de 10 h, colocación a 10 °C y manta de 40 mm. Se asignó peso uniforme a las ocho columnas de hidratación y se usaron los pesos configurados para GEFS frío, central y cálido. El límite evaluado es ΔT ≤ 20 °C y Tmax ≤ 70 °C; el criterio agregado de fallo es P_f ≤ 5 %.

- 20 de 24 combinaciones cumplen ambos límites de forma individual.
- Las 24 cumplen Tmax ≤ 70 °C; la mayor Tmax de la matriz es aproximadamente 56,21 °C.
- Cuatro exceden ΔT: M2_SS con clima frío (20,31 °C) y M4_SS con clima frío, central y cálido (21,45; 21,13; 20,80 °C).
- Ponderando los escenarios, P_T = 0 %, P_ΔT = P_f = 16,94 %. Por tanto, este plan no cumple el criterio agregado P_f ≤ 5 %.

La matriz incluye M1–M4 de las campañas SS y 22A. Como el IFC identifica la fundación como M1, M1_SS es el caso nominal más directamente alineado; las curvas M2–M4 deben interpretarse como sensibilidad de hidratación, no como mezclas confirmadas de esta fundación. M1_22A también es una curva de sensibilidad y tiene menor confianza indicada en el archivo fuente.

## Archivos de resultados

- [Matriz de 24 resultados por escenario](evidencia_semana_final/matriz_8_curvas_3_climas/resultados/resultados_por_escenario.csv)
- [Probabilidades ponderadas y cumplimiento agregado](evidencia_semana_final/matriz_8_curvas_3_climas/resultados/probabilidades_por_plan.csv)
- [Dossier con la configuración de la matriz](evidencia_semana_final/matriz_8_curvas_3_climas/dossier_corrida.json)
- [Serie térmica M1_SS / GEFS central con PUR 40 mm](evidencia_semana_final/manta_40mm/series_temperatura.csv)
- [Gráfico térmico del caso de referencia 40 mm](evidencia_semana_final/manta_40mm/evolucion_temperatura_40mm.png)
- [Ubicación de extremos del caso de referencia 40 mm](evidencia_semana_final/manta_40mm/ubicacion_extremos_40mm.png)
- [Comparación documentada de mantas de 30, 50 y 60 mm](evidencia_semana_final/documentacion/comparacion_mantas_30_50_60.md)
- [Informe del recálculo S25](evidencia_semana_final/documentacion/resumen_recalculo_s25.md)
- [Auditoría de salida E11 en UTF-8](evidencia_semana_final/auditoria_e11_utf8/resumen_e11.txt)

### Nota de ejecución

Las 24 simulaciones completaron y se guardaron en CSV. El evaluador dio un error al final al crear el manifiesto porque el dossier estaba fuera de la raíz del repositorio; por eso la carpeta de la matriz no contiene `manifest.json`. El error ocurrió después de calcular y exportar las métricas, no durante las corridas térmicas.
