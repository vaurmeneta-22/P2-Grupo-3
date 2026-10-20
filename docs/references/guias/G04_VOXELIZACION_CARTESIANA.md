# Guía G04 — Voxelización cartesiana

## Idea

No hacemos que la grilla siga la geometría.

Insertamos la geometría dentro de una grilla cartesiana:

```text
IFC/BIM
   ↓
superficies trianguladas
   ↓
clasificación inside/outside
   ↓
voxeles
```

Cada celda almacena como mínimo:

```python
active[i,j,k]
vertida_id[i,j,k]
material_id[i,j,k]
```

## Primera implementación

Para el baseline:

> La celda pertenece al sólido si el centro del voxel está dentro del sólido.

Eso produce una frontera staircase.

Es suficiente si se acompaña de un estudio de convergencia.

## Error de volumen

\[
V_h=N_{\mathrm{activo}}\Delta x\Delta y\Delta z
\]

\[
e_V=\frac{|V_h-V_{\mathrm{geom}}|}{V_{\mathrm{geom}}}
\]

Evaluar al menos tres resoluciones.

## Qué se debe aprender

La pregunta no es:

> “¿Se ve bonito?”

La pregunta es:

> “¿Cuánto cambia una magnitud relevante al refinar la discretización?”

Además del volumen, puede compararse posteriormente:

- \(T_{\max}\);
- \(\Delta T_{\max}\);
- tiempo al peak.

## Vertidas

Cada voxel debe conservar `vertida_id`.

La activación depende de:

```text
inicio_vertida[id]
T_colocacion[id]
```

Antes del inicio, el voxel no participa como hormigón activo.

## Fronteras y superficies

Una cara entre dos voxeles activos es una interfaz interna.

Una cara entre un voxel activo y el exterior puede recibir:

- convección;
- radiación;
- moldaje;
- aislación.

Una cara contra suelo puede conducir hacia material suelo modelado explícitamente o
utilizar una condición efectiva debidamente justificada.

## Nivel avanzado opcional

- fracciones de volumen;
- cut cells;
- signed-distance fields;
- octrees.

No son requeridos para aprobar.
