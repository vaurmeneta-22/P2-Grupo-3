# Guía G02 — Blender + Bonsai

## Rol de Blender

Blender es una suite libre y open source de creación 3D.

Referencia:
https://www.blender.org/about/

En este proyecto se utiliza principalmente como:

- viewport 3D;
- modelador geométrico;
- interfaz sobre la cual funciona Bonsai.

Blender por sí solo no convierte una malla en un objeto BIM.

## Rol de Bonsai

Bonsai agrega autoría BIM y trabajo IFC.

Instalación oficial:
https://docs.bonsaibim.org/quickstart/installation.html

La documentación actual indica instalar Blender y luego buscar **Bonsai** en:

```text
Edit → Preferences → Get Extensions
```

## Ejercicio mínimo

1. Crear/abrir un proyecto IFC.
2. Crear tres objetos distintos.
3. Asignar identidad BIM.
4. Nombrarlos:
   - `Vertida_01`
   - `Vertida_02`
   - `Vertida_03`
5. Guardar IFC.
6. Abrirlo desde Python con IfcOpenShell.
7. Modificar un objeto en Bonsai.
8. Guardar y comprobar que Python detecta el cambio.

## Property set sugerido

```text
Pset_HormigonMasivo

VertidaID
TipoHormigon
InicioPlanificado
TemperaturaColocacion_C
Moldaje
Aislacion
RetiroMoldaje_h
```

No es necesario que esos nombres sean idénticos si el grupo define y documenta un
schema equivalente.

## Qué no hacer

No usar el Python interno de Blender como núcleo del solver.

Queremos:

```text
Bonsai → IFC → Python normal
```

y no:

```text
Blender abierto → script interno → solver
```

## Recursos

Bonsai docs:
https://docs.bonsaibim.org/

Blender:
https://www.blender.org/

Video recomendado:
https://www.youtube.com/watch?v=MB68RUd98f8
