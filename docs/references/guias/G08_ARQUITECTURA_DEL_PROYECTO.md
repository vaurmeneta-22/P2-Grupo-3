# Guía G08 — Arquitectura del proyecto

## Pipeline

```text
Bonsai
  ↓
IFC
  ↓
IfcOpenShell
  ↓
geometry + properties
  ↓
voxelizer
  ↓
Cartesian grid
  ↓
thermal solver
  ↓
weather + hydration + construction
  ↓
Monte Carlo / HPC
  ↓
reliability
  ↓
construction plan
  ↓
VTK / PyVista / ParaView
```

## Estructura sugerida

```text
mass-concrete-planner/
├── README.md
├── pyproject.toml
├── geometry/
│   └── foundation.ifc
├── data/
│   ├── weather/
│   └── materials/
├── src/
│   └── mcplan/
│       ├── bim.py
│       ├── voxel.py
│       ├── thermal.py
│       ├── hydration.py
│       ├── boundary.py
│       ├── weather.py
│       ├── reliability.py
│       └── visualization.py
├── cases/
├── tests/
├── results/
└── docs/
```

## Principio 1 — Separar IO de física

Bueno:

```python
project = load_project(...)
grid = build_grid(project)
result = solve(grid, ...)
save_result(result)
```

Malo:

```python
# solver leyendo directamente botones/objetos de Blender
```

## Principio 2 — Inputs declarativos

El plan final debe poder modificarse sin cambiar el solver.

Ejemplo:

```yaml
vertidas:
  - id: 1
    start: "2026-11-27T04:00"
    placement_temperature_C: 16.0
    insulation_mm: 50
```

## Principio 3 — Tests antes de caso real

```text
analytical → verified solver → complex geometry
```

No al revés.

## Principio 4 — resultados reproducibles

Cada figura final debe tener:

- script o comando;
- input;
- commit;
- unidades;
- fecha del forecast cuando corresponda.

## Principio 5 — el modelo debe ser sometible a prueba

En la defensa el evaluador puede pedir cambiar:

- hora;
- temperatura de colocación;
- intervalo entre vertidas;
- aislación;
- retiro de moldaje.

El modelo debe responder sin reprogramación manual.
