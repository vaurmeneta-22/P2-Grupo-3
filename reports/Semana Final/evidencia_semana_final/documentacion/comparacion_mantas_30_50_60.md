# Comparación de mantas PUR

Se compararon espesores de 30, 50 y 60 mm, manteniendo iguales las demás condiciones: k = 0,025 W/(m·K), instalación en t = 15 h, retiro 360 h después de instalar, caras laterales y superior; colocación a 10 °C; mismo plan de vertidas, suelo, clima GEFS central, M1_SS, Δx = 0,25 m y horizonte de 336 h.

| Espesor | Tmax | Hora de Tmax | ΔT máximo | Hora de ΔT máximo | S25 frío en ese instante |
|---:|---:|---:|---:|---:|---:|
| 30 mm | 50,692 °C | 207,17 h | 19,749 °C | 333,5 h | 29,876 °C |
| 50 mm | 50,774 °C | 210,33 h | 16,730 °C | 336 h | 33,155 °C |
| 60 mm | 50,802 °C | 211,67 h | 15,832 °C | 336 h | 34,164 °C |

Con 30 mm, el ΔT queda **0,251 °C bajo el límite de 20 °C**. En el instante de su máximo, el voxel interior caliente está en `(-0,125; -0,125; 2,125) m`, a 49,625 °C; el punto S25 frío está en la cara lateral `x-`, en `(-8,975; -3,625; 1,125) m`, a 29,876 °C.

En los casos de 50 y 60 mm, el S25 más frío ocurre en la cara inferior hacia el suelo. S25 representa la temperatura reconstruida a 25 mm bajo la superficie, no el centro del voxel.

Los tres casos cumplen los límites de 70 °C y ΔT de 20 °C dentro de la ventana simulada. El máximo de 50 y 60 mm ocurre exactamente al final de las 336 h. El pronóstico también termina allí, así que esos casos no están comprobados después del día 14. El caso de 30 mm alcanza el máximo 2,5 h antes del final.

Cada carpeta `recalculo_30mm`, `recalculo_50mm` y `recalculo_60mm` contiene la configuración usada, manifiesto, CSV y temperaturas emparejadas en el máximo ΔT.

