# Sensibilidad al espesor de PUR del Plan 1

Se conservan las diez horas de colocación, el t0, los materiales, las caras de aplicación, la conductividad del PUR (0,025 W/m·K), su instalación a las 15 h, el suelo y los 24 escenarios ponderados del Plan 1. La manta se mantiene hasta las 672 h y el cálculo nominal usa malla de 0,25 m y paso máximo de 600 s. El único parámetro físico variado es el espesor del PUR.

El estudio prueba 80, 100, 120, 150, 200, 300 y 500 mm. Los 500 mm son un extremo hipotético de evaluación. El máximo instalable en obra dependería del producto y del detalle constructivo. El caso `pur_limite` usa un espesor equivalente de 10 m (`R = 400 m²·K/W`) para aproximar aislamiento casi perfecto en las caras cubiertas y visualizar el límite numérico. La condición del suelo en la cara inferior se conserva.

## Reproducción

```bash
for mm in 80 100 120 500; do
    python scripts/cierre_plan1.py "pur_${mm}mm" --procs 8
done
for variante in pur_150mm pur_200mm pur_300mm pur_limite; do
    python scripts/cierre_plan1.py "$variante" --procs 3 \
        --escenario M4_SS__GEFS_20261116_frio \
        --escenario M4_SS__GEFS_20261116_central \
        --escenario M4_SS__GEFS_20261116_calido
done
for mm in 100 120; do
    python scripts/cierre_plan1.py "pur_${mm}mm_fino" \
        --escenario M4_SS__GEFS_20261116_central
done
python scripts/resumir_cierre_plan1.py
```

El cálculo de probabilidad puntual requiere las 24 combinaciones. En variantes cribadas se informa exclusivamente la cota de la probabilidad que permiten los escenarios efectivamente calculados. La revisión en malla fina de 0,125 m se limita al escenario M4_SS central en cada espesor candidato.

## Resultados nominales

| Espesor | Escenarios | Diferencial máximo [°C] | Temperatura máxima [°C] | Riesgo ponderado a 672 h |
|---|---:|---:|---:|---:|
| 60 mm (referencia) | 24 | 21,503 | 56,298 | 25 % |
| 80 mm | 24 | 20,361 | 56,346 | 12,5 % |
| 100 mm | 24 | 19,645 | 56,378 | 0 % |
| 120 mm | 24 | 19,089 | 56,415 | 0 % |
| 150 mm | 3 | 18,524 | 56,455 | Cribado, 0/3 excedencias |
| 200 mm | 3 | 17,909 | 56,497 | Cribado, 0/3 excedencias |
| 300 mm | 3 | 17,287 | 56,563 | Cribado, 0/3 excedencias |
| 500 mm | 24 | 16,818 | 56,643 | 0 % |
| Aislamiento casi perfecto | 3 | 16,324 | 57,029 | Cribado, 0/3 excedencias |

Todos los resultados corresponden a la protección aún instalada. El diferencial se mantiene positivo aun en el límite numérico porque persisten la conducción dentro del hormigón y el intercambio térmico por la base. El procedimiento y la fecha de retiro para el nuevo espesor requieren una evaluación separada.

## Verificación y conclusión

Con el escenario M4_SS central y la malla refinada, el diferencial pasa de 19,645 a 19,928 °C con 100 mm, y de 19,089 a 19,372 °C con 120 mm. Ambos cumplen en esta comprobación, pero 120 mm deja 0,628 °C de margen frente a 0,072 °C con 100 mm. Las probabilidades de la tabla se refieren a la matriz nominal. Las 110 corridas nuevas tuvieron estabilidad y balance dentro de los límites comprobados; pasaron 211 pruebas automatizadas y hubo una omitida.

Los 500 mm alcanzaron su mayor diferencial exactamente a las 672 h. Ese valor es el máximo del intervalo simulado y no permite localizar el pico posterior. La reducción al pasar de 500 mm al límite numérico es de aproximadamente 0,49 °C en M4_SS central. El espesor de 120 mm es el candidato con mayor margen de los dos espesores verificados con malla fina para una evaluación posterior del retiro.
