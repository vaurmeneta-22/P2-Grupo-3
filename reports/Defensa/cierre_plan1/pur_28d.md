# Evaluación de PUR mantenido a 28 días

La decisión compara espesores de PUR en el Plan 1 de diez horas, con hormigón colocado a 10 °C e instalación de la manta a las 15 h. La manta permanece instalada durante las 672 h. La programación de las vertidas, el moldaje, el suelo, las tres trayectorias GEFS y las ocho curvas de hidratación se conservan. El paso máximo configurado es 600 s para ambas mallas. Se contabiliza una excedencia por escenario cuando la temperatura máxima supera 70 °C o el diferencial simultáneo con S25 supera 20 °C. La clasificación compara el peso exacto de la unión con 1/20.

Se reutilizan las cuatro matrices nominales de 60, 80, 100 y 120 mm en malla de 0,25 m. En la matriz fina de 100 mm se conserva la corrida de control M4_SS con clima central y se completan las otras 23 simulaciones de 672 h. Para contrastar 120 mm se reutilizan las 24 series refinadas que cubren 672 h. El script comprueba que las 673 muestras horarias del clima coinciden para cada miembro y que la serie refinada directa de M4_SS central coincide con el tramo recortado. Los extremos y las excedencias se recalculan sobre el corte a 672 h.

El resultado base de 336 h y los casos con retiro se conservan en sus archivos de origen. El retiro a las 375 h queda documentado en el [cierre del Plan 1](README.md).

## Reproducción

Desde la raíz del repositorio y con las dependencias del proyecto instaladas:

```bash
python scripts/cierre_plan1.py pur_100mm_fino --procs 8
python scripts/resumir_pur_28d.py
```

El ejecutor conserva un archivo comprimido por escenario y comprueba la firma de las entradas antes de reutilizarlo. Los archivos guardan los máximos con instantes y ubicaciones, los historiales S25, la estabilidad, el balance energético, el costo y la versión del entorno. En las series contrastadas de 120 mm, la estabilidad y el balance reportados corresponden a la corrida completa, y los máximos térmicos al corte de 672 h.

- [Resumen por espesor y malla](pur_28d_resumen.csv)
- [Comparación de 100 mm en los 24 escenarios](pur_28d_100mm_por_escenario.csv)
- [Comparación de 120 mm en los 24 escenarios](pur_28d_120mm_por_escenario.csv)
- [Figura de máximos S25](pur_28d_diferencial.png)

## Registro de ejecución

La campaña parte de `bd99654b27daae04dcd0b983fcb303223bdfd6e2`, que conserva la evidencia de `7266223989e6586d75abc5db5645780e0c45b224`. El piloto M4_SS cálido se ejecutó en un proceso y tomó 396 s. Los 22 escenarios restantes se ejecutaron con hasta ocho procesos. El entorno de cálculo fue Linux x86-64 en AMD EPYC 9V74, con nueve CPU lógicas disponibles, Python 3.12.14, NumPy 2.3.5 e IfcOpenShell 0.9.0.post1. La evaluación del repositorio con pytest registró 211 pruebas aprobadas y una omitida.

## Resultados y selección

| PUR instalado | Malla [m] | Temperatura máxima [°C] | Diferencial máximo S25 [°C] | Excedencias | Probabilidad conjunta |
|---|---:|---:|---:|---:|---:|
| 60 mm | 0,25 | 56,298 | 21,503 | 6/24 | 25 % |
| 80 mm | 0,25 | 56,346 | 20,361 | 3/24 | 12,5 % |
| 100 mm | 0,25 | 56,378 | 19,645 | 0/24 | 0 % |
| 100 mm | 0,125 | 56,437 | 19,928 | 0/24 | 0 % |
| 120 mm | 0,25 | 56,415 | 19,089 | 0/24 | 0 % |
| 120 mm | 0,125 | 56,470 | 19,372 | 0/24 | 0 % |

El PUR de 100 mm es el menor espesor evaluado que cumple en ambas mallas a 28 días. Con 100 mm y malla fina, M4_SS central alcanza el mayor diferencial a las 668 h 40 min y deja 0,072 °C de margen. La malla fina aumenta el máximo de la matriz en 0,283 °C. Con 120 mm el margen es de 0,628 °C. Como los máximos de ambas variantes refinadas ocurren en las últimas 24 h, la conclusión se circunscribe a las 672 h con PUR instalado. En la matriz fina de 100 mm, la mayor razón entre el paso aplicado y el límite estable fue 0,21334 y el mayor residuo relativo de balance fue inferior a 6×10⁻¹⁶.
