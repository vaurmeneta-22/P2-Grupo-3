# Cierre técnico del Plan 1 del Grupo 3

Continuación del árbol `7266223989e6586d75abc5db5645780e0c45b224`, verificado como HEAD de `main` al iniciar. Este identificador acredita la versión base de trabajo. El SHA presentado durante la defensa requiere evidencia de la entrega original.

## Entradas y criterio

El caso 4 se lee directamente de `configs/defensa/bateria7_casos.json`. Se mantienen la activación térmica instantánea de cada vertida en 0/2/4/7/9 h, las duraciones de 2/2/3/2/1 h, la colocación a 10 °C, el PUR de 60 mm en superficies expuestas desde 15 h hasta 375 h globales y el moldaje lateral de 18 mm en V1, V2 y V5 hasta las 72/74/81 h globales, respectivamente. Se conservan el suelo efectivo de 2 m a 20 °C, las propiedades térmicas y el IFC. La simulación llega a 672 h desde el 19-11-2026 a las 10:00 de America/Santiago (13:00 UTC).

Se reutiliza el GEFS descargado el 16-11-2026 en `data/weather/2026-11-16T213220Z_GEFS_ensemble/original.json`. La cola hasta el 17-12-2026 a las 10:00 local usa los mismos miembros completos seleccionados originalmente (07, control y 19), sin repetir ciclos ni empalmar miembros. Se comprobó que los instantes y temperaturas de las primeras 337 muestras horarias son idénticos para las tres trayectorias y que los tres metadatos extendidos identifican el SHA-256 del original. De 336 a 672 h, el pronóstico pierde habilidad y se interpreta como un escenario de estrés condicionado a la trayectoria archivada.

Cada una de las ocho curvas pesa 1/8, y los climas frío, central y cálido pesan 11/31, 10/31 y 10/31, respectivamente. Por independencia, cada combinación pesa 11/248 o 10/248. El indicador cuenta una sola vez cada escenario con T máxima estrictamente mayor que 70 °C o diferencial máximo estrictamente mayor que 20 °C. La suma se calcula como fracción exacta y se compara con 1/20. El diferencial resta, para cada instante, la mínima temperatura reconstruida a 25 mm de las caras del hormigón activo a la máxima temperatura simultánea del hormigón activo. En los eventos de retiro se registran los estados de borde anterior y posterior. El cambio en la reconstrucción S25 al variar la resistencia no implica un salto de temperatura en el interior del hormigón.

## Plazo constructivo y ajuste de retiro

La batería previa conserva sus once casos. Las alternativas protegidas de 10 y 24 h forman una comparación entre secuencias con duración física definida y buen desempeño hasta 336 h. Las 6,368 h del caso simultáneo representan la cota volumétrica a 100 m³/h. La activación conjunta de todas las vertidas en t = 0 es una idealización térmica. Un frente continuo necesita, en promedio, 63,68 m³/h para completar 636,776 m³ en diez horas. A 36 m³/h, el límite volumétrico es 17,69 h y excluiría también al Plan 1 de 10 h. La hipótesis de hasta 100 m³/h admite plazos inferiores. El repositorio señala que el caudal sostenible, número de bombas y continuidad de suministro siguen sin validación.

Se modeló una secuencia continua con un frente a 100 m³/h, duraciones obtenidas de los volúmenes IFC y llenado progresivo por niveles de voxel, con eventos cuantizados a 5 min mediante `tiempos_frente`. La temperatura de colocación, suelo, curvas y climas se mantienen. La manta sigue instalada de 15 a 375 h globales para aislar el cambio de cronograma y de activación. El moldaje conserva su regla de 72 h desde el inicio de su vertida. El cribado usa M4_SS con los climas frío y central. Sus cotas probabilísticas se derivan de los pesos completos, sin normalizar la muestra parcial.

Para estudiar el retiro se sustituyeron los 60 mm por tres capas equivalentes de 20 mm del mismo PUR, retiradas a 375, 480 y 576 h globales. La retirada gradual limita el salto inicial de resistencia y deja tiempo para reducir la diferencia núcleo–ambiente antes de exponer la superficie, pero también mantiene el núcleo caliente durante más tiempo. Esta configuración es un ensayo térmico, condicionado a la disponibilidad de un sistema de capas y a su aplicación en obra. Se cribó con los mismos dos escenarios críticos. La optimización de un protocolo de retiro basado en mediciones requiere validar protección, caudales, umbrales de intervención y condiciones climáticas en obra.

## Resultados y efecto sobre la decisión

| Ventana del Plan 1 exacto | Escenarios | T máxima [°C] | Diferencial máximo [°C] | Probabilidad conjunta |
|---|---:|---:|---:|---:|
| 0–336 h, batería publicada | 24 | 56,298 | 18,964 (336 h) | 0 % |
| 0–672 h, incluye retiro | 24 | 56,298 (214 h) | 36,112 (450 h 10 min) | 227/248 = 91,5323 % |

Se reprodujeron ambas métricas a tres decimales para los 24 escenarios de la ventana original, con discrepancia absoluta máxima de 0,000498 °C respecto de la tabla redondeada. En la extensión superaron 20 °C 22 escenarios y ninguno superó 70 °C. Los dos escenarios sin excedencia fueron M1_22A con clima frío y con clima cálido. Los máximos de T y del diferencial a 672 h corresponden a M4_SS con clima central. En el instante del diferencial máximo, el punto cálido está en (−0,125; −0,125; 2,125) m y el punto S25 más frío en (−8,975; −3,625; 1,125) m. Ninguno de los 24 escenarios tiene sus máximos de T o diferencial en las últimas 24 h del horizonte de 672 h. El cambio de clasificación se produce en la ventana posretiro y no depende de redondear la probabilidad.

La secuencia progresiva a 100 m³/h llegó a 35,195 °C (frío) y 36,089 °C (central) de diferencial máximo. El retiro escalonado llegó a 32,224 °C y 32,864 °C, respectivamente. Cada cribado determina una cota inferior de 21/248 = 8,4677 %, suficiente para descartar la configuración específica frente al 5 %. No se computa una probabilidad puntual a partir de solo dos escenarios. Los resultados no establecen que todo plazo inferior a diez horas sea inviable ni que otra regla de protección sea factible.

El Plan 1 conserva su resultado favorable a 14 días y deja de sostenerse como plan térmico hasta el retiro indicado. Para el informe posterior, la comparación debe distinguir los horizontes de 336 y 672 h, informar el retiro explícito y tratar los caudales como hipótesis. La Batería 7 y la propuesta B permanecen como evidencia independiente.

## Sensibilidad y comprobaciones

Se compararon dos refinamientos independientes en M4_SS con los climas frío y central, siempre hasta 672 h. Al reducir la malla a 0,125 m con 600 s, el diferencial subió 0,464–0,562 °C (frío, a 336 y 672 h) y 0,446–0,521 °C (central). La temperatura máxima aumentó cerca de 0,05 °C. Al reducir solo el paso máximo a 300 s con malla de 0,25 m, el diferencial cambió como máximo 0,0034 °C en valor absoluto a 672 h. Ambos escenarios siguen por encima del límite tras el refinamiento espacial y temporal, por lo que su peso conjunto de 21/248 por sí solo excede el 5 %. No se necesitó extrapolar el sesgo de malla a los otros escenarios para comprobar la clasificación. El 91,5323 % se refiere a la matriz nominal de 0,25 m y 600 s.

La peor razón paso aplicado/paso estable fue 0,213333, inferior al límite de estabilidad. El residuo energético relativo máximo se reporta en `resumen.csv`. Se ejecutaron 211 pruebas satisfactorias y una omitida, incluida una regresión que compara la simulación con y sin guardado de campos 3D.

## Reproducción y archivos

Desde la raíz del repositorio, con las dependencias del proyecto instaladas:

```bash
python scripts/cierre_plan1.py exacto --escenario M4_SS__GEFS_20261116_frio
python scripts/cierre_plan1.py exacto --procs 6
for variante in espacio tiempo progresivo_100 retiro_escalonado; do
  python scripts/cierre_plan1.py "$variante" --procs 2 \
    --escenario M4_SS__GEFS_20261116_frio \
    --escenario M4_SS__GEFS_20261116_central
done
python scripts/resumir_cierre_plan1.py
PYTHONPATH=src python -m pytest -q
```

Los comandos paralelos se ejecutaron en Linux. Por defecto, el lanzador trabaja en serie y admite volver a invocarse para completar solo los escenarios restantes. Cada registro `results/cierre_plan1/<variante>/<escenario>.json.gz` se escribe atómicamente y contiene la configuración efectiva, firma de entradas y código, fuentes de energía, eventos, extremos y serie por paso. Las tablas versionadas incluyen resultados por escenario con coordenadas, comparación a 336 h contra los valores publicados, sensibilidad, metadatos y tiempo por corrida. La figura y la tabla de curvas conservan el instante de los extremos; la figura traza el historial del escenario frío M4_SS.

Las corridas usan el mismo motor de diferencias finitas y mantienen los eventos y el muestreo horario de la batería original. La opción `guardar_campos=False` evita acumular matrices tridimensionales sin alterar integración, observación de S25, balance ni resúmenes. La estabilidad se verifica con la razón entre paso aplicado y límite del operador. Las comparaciones de espacio y tiempo cambian una variable a la vez. Se añaden las coordenadas del máximo interior concurrente con el máximo diferencial, además de la ubicación S25 y la del máximo global independiente.

Entorno de ejecución: Linux x86-64, CPU AMD EPYC 9V74, 9 CPU lógicas visibles y 9,7 GiB de memoria visibles, Python 3.12.14, NumPy 2.3.5, IfcOpenShell 0.9.0.post1 y Matplotlib 3.10.8. La prueba piloto nominal tardó 31,05 s y registró un máximo de 134 MiB de RSS. Los tiempos de matriz muestran tiempos de pared por proceso y su suma, que no equivale al tiempo transcurrido de la campaña en paralelo.

La continuación a 1008 h para 100 y 120 mm está documentada en [Evaluación posretiro](retiro_pur_100_120.md).
