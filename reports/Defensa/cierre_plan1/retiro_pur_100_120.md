# Extensión posretiro del Plan 1 con PUR de 100 y 120 mm

El corte previo de 672 h ya estaba evaluado con la manta instalada. Esta continuación añade el intervalo posretiro hasta 1008 h, conservando el caso original de 60 mm, los resultados previos de 100 y 120 mm y la batería de once casos. Mantiene el Plan 1 de 10 h, 10 °C de colocación y montaje de PUR a las 15 h.

## Meteorología y protocolos

El archivo original del GEFS del 16 de noviembre contiene valores comunes para los tres miembros seleccionados hasta la hora 733 desde el origen del plan. Se conservan los miembros 07, control y 19 completos hasta ese corte, con sus pesos originales. Las primeras 673 muestras son idénticas a las del corte de 28 días. Para las horas posteriores se usó el último ciclo de 24 h de cada miembro como continuación de estrés. Como contraste se mantuvo constante el último valor de cada miembro para el retiro en capas de 120 mm. La cola extrapolada es un supuesto de simulación y no un pronóstico.

Se comparó conservar la manta, el retiro total en la hora 672 y el retiro de cuatro capas iguales a las 672, 744, 816 y 888 h. El modelo suma resistencias de capas con la misma conductividad hasta su retiro. Las capas constituyen una configuración térmica hipotética cuya instalación requeriría verificación constructiva. El seguimiento cubre 336 h tras el retiro total o 120 h tras la última capa, hasta la hora 1008.

Los controles con manta mantenida se evalúan en los 24 escenarios nominales por espesor. Los retiros se criban en M4_SS frío y central, con pesos 11/248 y 10/248. Si ambos exceden el umbral, su peso combinado certifica que el protocolo incumple el criterio del 5 % en la matriz ponderada, aunque el riesgo total de ese protocolo no queda estimado.

## Reproducción

```bash
python scripts/preparar_clima_extendido_gefs.py --archivo-completo
for mm in 100 120; do
    python scripts/cierre_plan1.py "pur_${mm}mm_mantenido_1008" --procs 8
    for protocolo in total_672 capas_672; do
        python scripts/cierre_plan1.py "pur_${mm}mm_${protocolo}" --procs 2 \
            --escenario M4_SS__GEFS_20261116_frio \
            --escenario M4_SS__GEFS_20261116_central
    done
    python scripts/cierre_plan1.py "pur_${mm}mm_mantenido_1008_fino" \
        --escenario M4_SS__GEFS_20261116_central
done
python scripts/cierre_plan1.py pur_100mm_mantenido_1008_fino \
    --escenario M4_SS__GEFS_20261116_frio
for mm in 100 120; do
    python scripts/cierre_plan1.py "pur_${mm}mm_mantenido_1008_fino" \
        --escenario M4_SS__GEFS_20261116_calido
done
python scripts/cierre_plan1.py pur_120mm_capas_672_clima_constante --procs 2 \
    --escenario M4_SS__GEFS_20261116_frio \
    --escenario M4_SS__GEFS_20261116_central
python scripts/resumir_cierre_plan1.py
python scripts/graficar_retiro_pur.py
```

La carpeta de resultados contiene series y diagnósticos por escenario. Para preservar el alcance de las salidas anteriores, las variantes nuevas tienen nombres y carpetas distintos, y las originales conservan el horizonte de 672 h.

La búsqueda de duración conserva el cribado ya ejecutado para el frente progresivo de 6,368 h con el esquema de 60 mm. Esa comprobación excluye la configuración evaluada. La factibilidad de otras combinaciones más rápidas depende del caudal sostenible y de la protección elegida. La extensión de retiros mantiene fija la duración de diez horas.

## Resultados a 1008 h

| PUR | Protocolo | Escenarios | Máximo diferencial [°C] | Riesgo nominal o cota inferior |
|---|---|---:|---:|---:|
| 100 mm | Manta mantenida | 24 | 19,836 | 0 % |
| 120 mm | Manta mantenida | 24 | 19,496 | 0 % |
| 100 mm | Retiro total a 672 h | 2 | 32,404 | ≥ 8,4677 % |
| 120 mm | Retiro total a 672 h | 2 | 32,724 | ≥ 8,4677 % |
| 100 mm | Retiro en cuatro capas | 2 | 29,892 | ≥ 8,4677 % |
| 120 mm | Retiro en cuatro capas | 2 | 30,348 | ≥ 8,4677 % |
| 120 mm | Cuatro capas, cola climática constante | 2 | 29,393 | ≥ 8,4677 % |

El mayor diferencial del retiro escalonado de 120 mm con ciclo diario ocurre a 930 h 20 min, unas 42 h después de la última capa. La prueba con cola climática constante registra el máximo a 942 h 10 min. La conclusión de excedencia se mantiene para las dos continuaciones ensayadas. Los protocolos de retiro descritos quedan descartados por el criterio ponderado. El máximo registrado de temperatura se mantuvo bajo 70 °C en los casos ensayados.

## Sensibilidad y conclusión

Con la manta mantenida y la malla refinada de 0,125 m, el PUR de 100 mm alcanzó 20,077 °C con clima central, 20,101 °C con clima cálido y 19,964 °C con clima frío. El peso combinado de central y cálido es 20/248, igual a 8,0645 %. Es una cota inferior de riesgo en los escenarios refinados, suficiente para superar el criterio del 5 % en la hipótesis ensayada. Con 120 mm, las dos corridas refinadas llegaron hasta 19,739 °C, con margen mínimo de 0,261 °C en los escenarios verificados. Se mantiene como candidato con manta puesta para esta ventana. El retiro total y escalonado probados quedan descartados. Continúa abierta la selección de un protocolo de retiro alternativo.

Las 63 corridas nuevas conservaron el balance relativo por debajo de 1,1×10⁻¹⁵ y el paso aplicado por debajo del límite estable. Se comprobó igualdad exacta del historial nominal en todas las muestras anteriores a 672 h respecto de la corrida correspondiente de 28 días. Se ejecutaron 211 pruebas automatizadas y hubo una omitida. Los tiempos de corrida y detalles de equipo están en los registros de ejecución.
