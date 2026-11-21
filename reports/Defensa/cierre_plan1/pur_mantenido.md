# Plan 1 con PUR mantenido durante 28 días

Se extendió el Plan 1 de diez horas con la manta PUR de 60 mm instalada a las 15 h y sin retiro dentro de la simulación. Se conservaron la programación de las vertidas, las propiedades del hormigón, las condiciones de borde restantes, las ocho curvas, el GEFS extendido con los mismos tres miembros y la distribución de pesos. El horizonte de 672 h es común con la variante original, que retira la manta a las 375 h. Se usa la malla nominal de 0,25 m y el paso máximo de 600 s.

```bash
python scripts/cierre_plan1.py pur_mantenido --procs 8
python scripts/resumir_cierre_plan1.py
```

Los registros por escenario se guardan en `results/cierre_plan1/pur_mantenido/`. Las tablas distinguen el escenario con manta mantenida del retiro original y de las otras variantes. Las celdas de diagnóstico de retiro quedan vacías cuando el PUR permanece puesto.

## Resultados

| Condición de la manta | Horizonte | Máximo S25 | Probabilidad conjunta ponderada |
|---|---:|---:|---:|
| Mantenida | 336 h | 18,964 °C | 0 % |
| Retirada a las 375 h | 672 h | 36,112 °C a las 450 h 10 min | 91,5323 % |
| Mantenida | 672 h | 21,503 °C a las 597 h 20 min | 25 % |

Las series con y sin retiro coinciden exactamente antes de las 375 h en los 24 escenarios. Con la manta mantenida, el primer cruce individual de 20 °C sucede a las 397 h (peso de 11/248, 4,4355 %). A las 402 h 30 min se incorpora un segundo escenario y el acumulado alcanza 21/248 (8,4677 %). A las 672 h, seis escenarios han superado 20 °C y ninguno supera 70 °C. La temperatura máxima es 56,298 °C. La extensión reduce sustancialmente la excedencia en comparación con el retiro original, pero todavía excede el criterio de 5 %. El protocolo de retiro posterior a 672 h queda fuera de la evaluación.
