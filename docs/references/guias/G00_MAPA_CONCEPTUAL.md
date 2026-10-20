# Guía G00 — Mapa conceptual

## BIM no es un programa

**BIM (Building Information Modeling)** es una manera estructurada de representar y
gestionar información del activo construido.

Puede incluir:

- identidad;
- geometría;
- materiales;
- propiedades;
- relaciones espaciales;
- fases;
- cantidades;
- información constructiva.

Un modelo BIM útil para este proyecto puede expresar:

```text
Fundación
├── Vertida_01
│   ├── geometría
│   ├── hormigón = M1
│   ├── T_colocación = 18 °C
│   ├── inicio = ...
│   └── aislación = ...
├── Vertida_02
└── ...
```

## IFC

**IFC (Industry Foundation Classes)** es un estándar abierto e internacional de datos
para describir el entorno construido.

buildingSMART lo define como una descripción digital estandarizada y vendor-neutral.

Referencia oficial:

https://technical.buildingsmart.org/standards/IFC/

IFC puede contener:

- objetos;
- GUIDs;
- tipos;
- geometría;
- materiales;
- property sets;
- relaciones;
- ubicación.

IFC **no** es:

- un solver;
- un motor de diferencias finitas;
- una API meteorológica;
- un formato óptimo para millones de valores de temperatura.

## Capas del proyecto

```text
AUTORÍA BIM
Bonsai / Blender
      ↓
IFC
      ↓
IfcOpenShell
      ↓
GEOMETRÍA COMPUTACIONAL
voxelización cartesiana
      ↓
FÍSICA
solver térmico propio
      ↓
INCERTIDUMBRE
Monte Carlo / HPC
      ↓
RESULTADOS
VTK
      ↓
VISUALIZACIÓN
PyVista / ParaView
      ↓
DECISIÓN
plan de vertidas
```

La separación conceptual más importante del proyecto es:

> **BIM ≠ solver ≠ visualizador.**
