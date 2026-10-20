# Guía G01 — BIM, IFC, Revit y openBIM

## Revit

Autodesk Revit es una aplicación de autoría BIM ampliamente utilizada en la industria.

Su formato nativo es propietario (`.rvt`).

Revit puede importar/exportar IFC, permitiendo un flujo:

```text
Revit → IFC → software independiente
```

Nuestro proyecto no dependerá de Revit.

## Bonsai

Bonsai permite trabajar con BIM/IFC desde Blender utilizando un stack abierto.

```text
Blender = interfaz y modelado 3D
Bonsai = funcionalidad BIM/IFC
IfcOpenShell = motor/biblioteca IFC
```

## openBIM

La idea de openBIM es evitar que el intercambio dependa de un único fabricante.

En este proyecto el contrato será:

```text
herramienta de autoría
        ↓
       IFC
        ↓
nuestro código
```

Por ello, conceptualmente, un IFC bien exportado desde Revit, Bonsai u otra herramienta
debería poder alimentar el mismo pipeline.

## Interoperabilidad no significa perfección

Un IFC debe validarse:

- unidades;
- clases;
- geometría;
- propiedades;
- placements;
- nombres/GUIDs;
- property sets.

Dos programas distintos pueden exportar el mismo concepto de maneras diferentes.

## Referencias

buildingSMART — IFC:
https://technical.buildingsmart.org/standards/IFC/

Bonsai:
https://docs.bonsaibim.org/

IfcOpenShell:
https://docs.ifcopenshell.org/

Autodesk — BIM interoperability:
https://www.autodesk.com/industry/bim/interoperability

## Video inicial

Getting Started with Bonsai BIM — IFC Addict:
https://www.youtube.com/watch?v=MB68RUd98f8
