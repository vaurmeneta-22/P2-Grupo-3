# E10 · Entorno y clima

Qué agrega E10 al pipeline IFC 3D y cuáles son sus supuestos. Todo está **desactivado por
defecto** (`entorno.activo = false`): con la configuración anterior el pipeline reproduce E09
byte a byte (control térmico, historial y métricas E09).

## 1. Condiciones de borde por tipo de cara (`mcplan/thermal/bordes.py`)

| Cara | Modelo |
|---|---|
| Expuesta al aire | Robin: `q = (T_aire(t) − T_cara) / R_total`, `h` configurable (`entorno.aire.h_w_m2_k`) |
| Base física (cara `z−` en el plano más bajo del hormigón) | Suelo efectivo (ver abajo) |
| Con moldaje o aislación | Resistencias en serie con media celda y convección |
| Entre dos celdas activas | Sin cambios: conducción por conductancia (media armónica de k) |

Un voladizo con cara `z−` expuesta que no esté en la base usa aire.

### Suelo efectivo
Capa de espesor `L` y conductividad `k_s` hacia una temperatura profunda `T_suelo`:
`h_eq = k_s / L`, equivalente a un Robin con `T_amb = T_suelo`. Parámetros en
`entorno.suelo`: `temperatura_c`, `conductividad_w_m_k`, `espesor_efectivo_m`
(valores de la demo: 18 °C, 1,5 W/(m K), 0,8 m; **supuestos, no medidos**).

Supuestos y límites:
- Régimen casi estacionario en el suelo: no hay acumulación de calor en él. Al inicio
  sobreestima la pérdida (el suelo real se calienta) y la profundidad de penetración en 14 días
  es del orden de `√(α t)` ≈ 0,8 m (α ≈ 5·10⁻⁷ m²/s), por eso `L = 0,8 m` por defecto.
- Temperatura profunda constante y uniforme; sin humedad, napa ni capa de sello.
- Flujo positivo entra al hormigón: con suelo más frío que la base el flujo es negativo.

### Moldaje y aislación como resistencias en serie
```
R_total = d/(2 k_hormigón) + Σ R''_capas + 1/h        flujo = ΔT / R_total   [W/m²]
```
`R'' = espesor / conductividad` (o `resistencia_m2_k_w` directa). Varias capas sobre la misma
cara se suman. La temperatura de cara estimada (seguimiento y extremos) sigue siendo la del
**hormigón**, es decir, la cara bajo la protección.

### Protecciones y retiro, por vertida (`entorno.protecciones_por_vertida`)
```jsonc
"1": [ { "nombre": "moldaje", "caras": ["x-","x+","y-","y+"],
         "espesor_m": 0.018, "conductividad_w_m_k": 0.13,
         "retiro_h": 72.0, "retiro_desde": "colocacion" } ]
```
- `retiro_desde`: **`colocacion`** (por defecto) mide desde `inicio_h` de *esa* vertida;
  **`inicio_obra`** mide desde t = 0 de la simulación. Un retiro anterior a la colocación se rechaza.
- `retiro_h`: número, `null` (no se retira) o `"ifc"` (usa `RetiroMoldaje_h` del IFC, 72 h; el IFC no
  documenta su referencia temporal, por eso se interpreta con `retiro_desde`).
- Las propiedades IFC `Moldaje` y `Aislacion` están en `TBD`: los espesores y conductividades
  de la demo son valores genéricos de referencia, no del proyecto.
- El retiro es un instante de evento en `simular_vertidas`: se integra hasta él, se actualizan las
  conductancias **sin tocar la temperatura** y el intervalo siguiente recalcula el paso estable
  (`ResultadoEtapas.retiros` guarda el paso antes y después). Con la malla de referencia el
  paso de 600 s es el que manda (el estable es ~45 000 s), así que el recálculo no cambia el paso.

## 2. Clima (`mcplan/weather/`)

- `loader.construir_aire(config["clima"], ...)` entrega `T_aire(t)` para el Robin
  (`t` en segundos desde `inicio_obra_utc`). Fuentes: `constante`, `sintetica`, `snapshot`.
- **La serie de la demo es sintética** y está identificada como tal (`aire.es_sintetica`,
  `entorno_usado.json`). No es un pronóstico.
- Sitio del Grupo 3: **Ovalle / Ecos del Mar, Coquimbo**. Latitud y longitud están **pendientes de
  `sitios.yaml`** (`sitio.latitud/longitud = null`). Eso solo bloquea la consulta real:
  `descargar_y_guardar` lanza `SitioSinCoordenadas` antes de tocar la red.
- Primer cliente: Open-Meteo (`openmeteo.py`), `temperature_2m` horaria en UTC. URL y parseo se
  prueban sin red.

### Viento y radiación sin contar dos veces el mismo intercambio (no implementado)
- **Viento**: se modela con un único `h(v) = h_natural + b·v` combinado (convectivo + radiativo
  linealizado) que reemplaza al `h` constante; no se agrega un segundo coeficiente.
- **Radiación solar**: entra como flujo absorbido `q_sol = α·G` [W/m²] (`α` = absorbancia de la
  superficie, `G` = radiación sobre ella), con `T_aire` en el Robin; o, de forma equivalente, como
  temperatura aire-sol `T_as = T_aire + α·G / h` en el *mismo* Robin. Nunca ambas, y no se
  suma si el intercambio ya es un efectivo (aire-sol o `h` calibrado con sol incluido).
- **Onda larga**: si `h` ya es combinado, no se suma un término de onda larga aparte.

## 3. Snapshots de pronóstico (`mcplan/weather/snapshots.py`)

Cada descarga crea una carpeta nueva `data/weather/<descarga UTC>/` (si el nombre existe se
agrega `_2`, `_3`; no se sobrescribe ni se borra nada):

| Archivo | Contenido |
|---|---|
| `original.json` | bytes exactos recibidos |
| `forecast.json` | serie normalizada: `time_utc` (ISO, UTC), `temperature_2m` en degC |
| `metadata.json` | `source`, `model`, `site`, `latitude`, `longitude`, `timezone`, `issued_at` (null si la fuente no lo informa), `retrieved_at`, `valid_from`, `valid_to`, `units`, `original_file`, `original_sha256`, `synthetic` |

`metadata.json` incluye `series_kind`: `pronostico`, `reanalisis`, `observacion` o `sintetica`.
`leer_snapshot` verifica unidades y el sha256 del original. `normalizar_serie` rechaza tiempos
duplicados o desordenados, unidades desconocidas (acepta degC y K) y datos faltantes (o los
descarta si se pide, y la cobertura real se informa aparte).

**Cobertura:** `cobertura(serie, inicio, fin)` informa el intervalo disponible y las horas
faltantes. Si el pronóstico no cubre el horizonte, `construir_aire` falla y muestra ese intervalo,
salvo que se pida `snapshot.extension` (`ciclo_diario` repite el último ciclo de 24 h;
`ultimo_valor` mantiene el último dato). Cada punto agregado queda marcado como
`extension:<método>` en la serie y en `descripcion["puntos_por_origen"]`: es una hipótesis del
equipo, no un pronóstico.

## 4. Cómo ejecutar

- Demo corta: abrir `scripts/demo_e10_entorno.py` y ▶️ (configuración `configs/e10_demo.jsonc`;
  salida `results/e10_demo/demo_e10.png`).
- Pipeline completo con el entorno: copiar los bloques `sitio`, `clima` y `entorno` de
  `configs/e10_demo.jsonc` a la configuración deseada.
- Pruebas: `python -m unittest discover -s tests -p "test_*.py"` con `PYTHONPATH=src`.
