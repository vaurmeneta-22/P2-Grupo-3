# E11 · Clima desde snapshot y exportación 3D

Reutiliza lo de E10 (bordes, loader, snapshots). Estado: la demo usa un **pronóstico real de Open-Meteo**
(`best_match`) para el punto del Grupo 3 (Mantos de Hornillos / Ovalle, Camarico – Ecos del Mar, Coquimbo):
latitud -31.1115, longitud -71.6147, tomadas de la tabla docente de sitios (procedencia en
`coordinates_source` de los metadatos; `sitios.yaml` no estaba disponible para contrastarlas).
El proveedor devolvió el punto de grilla -31.107206, -71.615204 a 326 m de elevación (queda en `original.json`).
El snapshot sintético sigue disponible en `data/weather/2026-11-10T160000Z/` (apuntar `clima.snapshot.carpeta`
ahí para ensayar sin red).

## Archivos

| Qué | Dónde |
|---|---|
| Config E11 (se superpone a la base y a E10) | `configs/e11_clima_real.jsonc` |
| Demo (▶️) | `scripts/demo_e11_clima_real.py` |
| Visor 3D (PyVista) | `scripts/ver_resultado_3d.py` (instalar con `requirements-visor.txt`) |
| Snapshot real (pronóstico, descargado 2026-11-10 16:12 UTC) | `data/weather/2026-11-10T161233Z/` |
| Snapshot sintético | `scripts/crear_snapshot_sintetico_e11.py` → `data/weather/2026-11-10T160000Z/` |
| Descarga | `scripts/descargar_snapshot_e11.py` |
| Exportador VTK | `src/mcplan/visualization/vtk_export.py` |
| Pruebas | `tests/test_e11_vtk_clima.py` |

## Cómo repetir o renovar la descarga

1. Las coordenadas están en el bloque `sitio` de `configs/e11_clima_real.jsonc` (o se sustituyen con `--lat/--lon`).
2. `python scripts/descargar_snapshot_e11.py --modo pronostico` (hasta 16 días) o
   `--modo historico --inicio AAAA-MM-DD --fin AAAA-MM-DD` (reanálisis).
3. Copiar la carpeta impresa a `clima.snapshot.carpeta`. Cada descarga crea una carpeta nueva con
   `original.json`, `forecast.json` y `metadata.json` (proveedor, producto, coordenadas, unidades, descarga,
   intervalo cubierto, `series_kind`; `issued_at` queda `null` porque Open-Meteo no informa la emisión).

## Dónde se lee el aire (T_aire(t))

`loader.construir_aire` → `SerieTemperatura.a_funcion` (interpolación lineal; `t` en s desde el inicio de obra)
→ `robin(h, aire.funcion)` en `crear_bordes_entorno` → `DifusionCartesiana._flujo_borde_w`
(`cond.temperatura(tiempo_s)`). La demo imprime archivo y línea de cada paso.
Ecuación por cara expuesta: `Φ = G · (T_aire(t) − T_celda)`, `G = A / (d/(2k) + R''_protección + 1/h)`.
`h` es constante y configurable (`entorno.aire.h_w_m2_k`). Viento y radiación **no están activados**;
su incorporación (un único `h(v)` combinado y radiación como flujo absorbido `α·G` o temperatura aire-sol,
nunca ambos) está descrita en `docs/E10_entorno_y_clima.md`.

## Formato 3D

`results/e11_demo/vtk/estado.pvd` (colección temporal, `timestep` en **horas**) y un `.vtu` por estado.
Solo se escriben los voxeles **activos**. Coordenadas en **metros**. Datos de celda: `temperatura_C` (°C),
`vertida_id`, `mascara_activa`, `n_caras_expuestas`. Datos de campo: `TimeValue`/`tiempo_h` (h), `tiempo_s` (s).
`marcas.pvd` + `marcas_*.vtu`: puntos con `codigo_marca` = 1 núcleo, 2 voxel próximo a la cara (T del voxel),
3 cara exterior con **T reconstruida** (resistencia de media celda y borde del solver; no es medición).

## Guion del lab (jueves)

1. Mostrar `data/weather/2026-11-10T161233Z/metadata.json` (y `original.json`, `forecast.json`).
2. Abrir `results/e11_demo/aire_snapshot.png` (la demo lo regenera).
3. Abrir `scripts/demo_e11_clima_real.py` y ▶️ (reiniciar la ventana si ▶️ falla: ver `.vscode/settings.json`).
4. Abrir el resultado: `python scripts/ver_resultado_3d.py` (PyVista) o ParaView con `estado.pvd` + `marcas.pvd`.
5. Cortar: mover el plano de corte (o `--clip`); señalar el núcleo (esfera roja) y la cara exterior reconstruida (azul).
   `--captura DIR` genera `3d_completo.png`, `3d_slice.png`, `3d_clip.png` sin ventana.

## Ventana de la corrida real

`inicio_obra_utc = 2026-11-11T12:00:00Z`, 96 h (hasta 2026-11-15 12:00Z), dentro de la cobertura del snapshot
(2026-11-10 00:00Z a 2026-11-25 23:00Z, 384 puntos horarios, 12,7 a 21,4 °C). Es un **pronóstico**: no es una observación.
Open-Meteo no informa la hora de emisión (`issued_at = null`).
