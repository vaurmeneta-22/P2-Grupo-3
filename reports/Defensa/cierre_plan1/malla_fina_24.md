# Matrices de malla fina para la recomendación del Plan 1

Se comparan el Plan 1 original (60 mm, 336 h) y la alternativa extendida (120 mm con la protección mantenida, 1008 h). Ambas matrices utilizan ocho curvas, tres historias climáticas y malla de 0,125 m con paso máximo de 600 s. La tabla consolidada calcula las probabilidades con los pesos completos. Los casos finos ya ejecutados de 120 mm central y cálido se conservan.

```bash
python scripts/cierre_plan1.py plan1_336_fino --procs 8
python scripts/cierre_plan1.py pur_120mm_mantenido_1008_fino --procs 8
python scripts/resumir_malla_fina_plan1.py
```

La comprobación de retiro previa evaluó retiro total a 672 h y retiro en cuatro capas a 672, 744, 816 y 888 h. Ambos protocolos excedieron el criterio con 100 y 120 mm en el cribado. Cada par de escenarios con excedencia pesa 21/248. Una secuencia distinta sería otro protocolo, por ejemplo variar el tamaño de cada capa y el intervalo entre retiros. El retiro se conserva como observación técnica y no se usa para calificar la ventana con manta instalada.

| Configuración | Escenarios con excedencia | Temperatura máxima | Diferencial máximo | Riesgo ponderado |
|---|---:|---:|---:|---:|
| 60 mm, 336 h | 0 de 24 | 56,347 °C | 19,428 °C | 0 % |
| 120 mm mantenido, 1008 h | 0 de 24 | 56,470 °C | 19,739 °C | 0 % |

El Plan 1 original conserva 0,572 °C de margen en el diferencial y la alternativa extendida conserva 0,261 °C, bajo las hipótesis climáticas modeladas. Se ejecutaron 46 corridas nuevas y se reutilizaron dos de la misma configuración. La razón máxima de estabilidad fue 0,21334 y el mayor residuo relativo del balance fue inferior a 6,5×10⁻¹⁶. Los datos comprimidos conservan los máximos y sus ubicaciones por escenario.
