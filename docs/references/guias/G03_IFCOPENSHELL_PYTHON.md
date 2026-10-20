# Guía G03 — IfcOpenShell desde Python

## Objetivo

Transformar:

```text
IFC
```

en:

```text
objetos Python + geometría triangulada + propiedades
```

## Abrir un archivo

```python
import ifcopenshell

model = ifcopenshell.open("modelo.ifc")
print(model.schema)
```

Tutorial oficial:
https://docs.ifcopenshell.org/ifcopenshell-python/hello_world.html

## Consultar elementos

```python
for obj in model.by_type("IfcProduct"):
    print(obj.is_a(), obj.Name, obj.GlobalId)
```

El `GlobalId` es útil como identidad persistente.

## No parsear geometría IFC manualmente

IFC permite representar geometría mediante:

- extrusiones;
- perfiles;
- booleanas;
- BReps;
- mallas;
- placements anidados.

IfcOpenShell ya resuelve estas representaciones.

Referencia:
https://docs.ifcopenshell.org/ifcopenshell-python/geometry_processing.html

## Geometry iterator

Para procesar muchos elementos, la documentación recomienda el geometry iterator:

https://docs.ifcopenshell.org/ifcopenshell/geometry_iterator.html

Este puede entregar:

```text
vertices
faces trianguladas
transformación
materiales
GUID
```

y soporta procesamiento multicore/caching.

## Capa de abstracción recomendada

El módulo BIM debería devolver algo simple:

```python
VertidaGeometry(
    guid=...,
    name="Vertida_01",
    vertices=...,
    faces=...,
    properties=...,
)
```

Luego `voxel.py` no necesita conocer detalles internos del esquema IFC.

## Test mínimo

```bash
python inspect_ifc.py modelo.ifc
```

Debe reportar:

```text
Vertida       GUID       Inicio       Tcolocación
Vertida_01    ...        ...          ...
...
```

## Recursos

IfcOpenShell:
https://docs.ifcopenshell.org/

Code examples:
https://docs.ifcopenshell.org/ifcopenshell-python/code_examples.html

Geometry processing:
https://docs.ifcopenshell.org/ifcopenshell-python/geometry_processing.html

Webinar relacionado con Python + IfcOpenShell:
https://www.youtube.com/watch?v=j2GMbV9vkm8
