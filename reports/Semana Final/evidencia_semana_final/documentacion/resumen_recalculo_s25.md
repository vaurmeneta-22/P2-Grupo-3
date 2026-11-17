# Recalculo de S25 — comparación 50 mm / 60 mm

## Corrección

El observador S25 guardaba en caché las superficies expuestas. Al activarse una vertida o cambiar una protección, podía seguir midiendo una cara que ya era una interfaz interior. Se corrigió la firma de caché para que incluya la máscara activa, las resistencias de borde y las selecciones de superficie por cara.

## Caso comparado

- Inicio: 19-11-2026 10:00, hora de Ovalle.
- Horizonte: 336 h; `dx = 0.25 m`, `dt máximo = 600 s`.
- Plan: V1 0 h/2 h; V2 2 h/2 h; V3 4 h/3 h; V4 7 h/2 h; V5 9 h/1 h.
- Temperatura de colocación: 10 °C en las cinco vertidas.
- Fuente: M1_SS; clima GEFS central.
- Suelo: capa efectiva de 2 m, 20 °C, `k = 1.5 W/(m K)`.
- Moldaje: V1, V2 y V5, `k = 0.1 W/(m K)`, 18 mm.
- Manta PUR: `k = 0.025 W/(m K)`, instalada en `t = 15 h`, retiro 360 h después de instalar. Solo se cambió el espesor entre los casos.

## Resultados recalculados

| Manta | Tmax (°C) | Hora de Tmax | ΔT máximo (°C) | Hora de ΔT máximo | S25 frío (°C) en ese instante |
|---|---:|---:|---:|---:|---:|
| 50 mm | 50.774 | 210.33 h | 16.730 | 336 h | 33.155 |
| 60 mm | 50.802 | 211.67 h | 15.832 | 336 h | 34.164 |

En ambos casos el punto S25 más frío está en la cara inferior `z-`, en `(-8.875, -3.625, 0.025) m`. Es una temperatura reconstruida 25 mm bajo la superficie, no el centro de una celda. El centro de la celda de borde más cercana está en `(-8.875, -3.625, 0.125) m`.

A `t = 336 h`, la pareja que define ΔT para la manta de 60 mm es: voxel interior caliente `49.997 °C` en `(-0.125, -0.125, 2.375) m`, menos S25 frío `34.164 °C` en la base = `15.832 °C`.

## Comprobaciones y límite

- Suite automatizada: **184 pruebas aprobadas**.
- El test nuevo reproduce el cambio de S25 al activarse otra vertida y al retirar una protección.
- Los dos casos cumplen `ΔT ≤ 20 °C` y `Tmax ≤ 70 °C` dentro de los 14 días simulados.
- El máximo ΔT cae exactamente en el último instante del horizonte. El snapshot climático termina en esas mismas 336 h; por tanto, no se ha comprobado si el ΔT permanece bajo 20 °C después del día 14.

Los JSON reproducibles, CSV, manifiestos y la pareja de temperaturas están en sus carpetas `recalculo_50mm` y `recalculo_60mm`. Los JSON guardados deben copiarse a `configs/` dentro del repositorio para volver a ejecutar, porque sus rutas de entradas son relativas a la raíz del proyecto.

