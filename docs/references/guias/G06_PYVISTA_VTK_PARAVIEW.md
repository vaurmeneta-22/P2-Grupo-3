# Guía G06 — PyVista, VTK y ParaView

## Por qué no Matplotlib para todo

Matplotlib sigue siendo ideal para:

- curvas;
- convergencia;
- histogramas;
- barras;
- sensibilidad.

Pero para campos 3D necesitamos:

- slices;
- clipping;
- isosuperficies;
- transparencias;
- animación;
- grandes datasets.

## VTK

VTK es el ecosistema/modelo de datos científico subyacente.

No es necesario programar directamente contra la API de VTK.

## PyVista

PyVista ofrece una API Python mucho más simple sobre VTK.

Getting Started:
https://docs.pyvista.org/getting-started/

Tutorial:
https://tutorial.pyvista.org/

Para una grilla uniforme, `pyvista.ImageData` representa naturalmente dimensiones,
spacing y origin.

Documentación:
https://docs.pyvista.org/api/core/_autosummary/pyvista.imagedata

Ejemplo conceptual:

```python
import pyvista as pv

grid = pv.ImageData(
    dimensions=(nx+1, ny+1, nz+1),
    spacing=(dx, dy, dz),
    origin=(xmin, ymin, zmin),
)

grid.cell_data["temperature"] = T.ravel(order="F")
grid.cell_data["vertida_id"] = vertida_id.ravel(order="F")

grid.save("temperature_0048h.vti")
```

PyVista soporta guardar `ImageData` como `.vti`.

## ParaView

ParaView es una aplicación open source de visualización científica.

Documentación:
https://docs.paraview.org/en/latest/

Tutorial:
https://docs.paraview.org/en/latest/Tutorials/SelfDirectedTutorial/introduction.html

Puede manejar series temporales de archivos y grandes datasets.

## Outputs recomendados

```text
temperature_0000.vti
temperature_0001.vti
...
```

o una serie/metaarchivo VTK compatible.

## Visualizaciones mínimas

- temperatura 3D;
- slice vertical;
- slice horizontal;
- trayectoria núcleo-superficie;
- vertida_id;
- instante y ubicación de \(T_{\max}\);
- mapa de excedencia;
- probabilidad espacial de excedencia en la etapa de confiabilidad.

## Video recomendado

ParaView Webinar: First Steps — Kitware:
https://www.youtube.com/watch?v=IjXnZ8WeEJw
