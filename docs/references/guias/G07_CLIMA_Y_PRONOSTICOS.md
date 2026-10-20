# Guía G07 — Clima, pronósticos y trazabilidad

## Fuente recomendada

Open-Meteo:

https://open-meteo.com/

APIs útiles:

- Forecast:
  https://open-meteo.com/en/docs
- Ensemble:
  https://open-meteo.com/en/docs/ensemble-api
- Historical Forecast:
  https://open-meteo.com/en/docs/historical-forecast-api
- Previous Runs:
  https://open-meteo.com/en/docs/previous-runs-api
- Historical Weather:
  https://open-meteo.com/en/docs/historical-weather-api

## Variables

### Obligatoria

- temperatura del aire a 2 m.

### Si se modela convección dependiente del viento

- velocidad del viento a 10 m.

### Si se modela radiación solar

- shortwave radiation / GHI;
- eventualmente direct/diffuse según el modelo.

No agregar humedad simplemente porque esté disponible: debe existir una ecuación donde
esa variable participe.

## Forecast no es observación

Debe distinguirse entre:

```text
observación / reanalysis
forecast determinista
ensemble forecast
```

La fecha de descarga importa.

Un pronóstico descargado una semana antes no es el mismo dato que el pronóstico
descargado el día anterior.

## Snapshots

Ejemplo:

```text
data/weather/
├── 2026-11-11T1200/
│   ├── metadata.json
│   └── forecast.json
├── 2026-11-20T1200/
└── FINAL_TBD/
```

Metadata mínima:

```json
{
  "source": "...",
  "model": "...",
  "retrieved_at": "...",
  "latitude": "...",
  "longitude": "...",
  "timezone": "...",
  "variables": ["temperature_2m", "..."],
  "units": {"temperature_2m": "degC"}
}
```

## Horizonte incompleto

Antes de la fecha final (TBD) el forecast puede no cubrir los 14 días posteriores a la
construcción.

Eso es parte del problema, no un error del curso.

El equipo debe definir cómo completa la cola faltante:

- clima histórico;
- distribución condicional;
- ensemble de mayor alcance;
- escenario conservador;
- otro método.

## Ensemble y confiabilidad

Cuando existan miembros de ensemble útiles:

\[
W^{(1)}(t),W^{(2)}(t),...,W^{(N)}(t)
\]

pueden entrar directamente como escenarios meteorológicos.

No es obligatorio convertir cada hora del forecast en una variable normal independiente.

## Validación útil opcional

La API de previous runs/historical forecast permite analizar cómo cambió el forecast
según el lead time y puede ayudar a caracterizar error meteorológico.
