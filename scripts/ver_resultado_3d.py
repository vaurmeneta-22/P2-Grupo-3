"""Visor 3D (PyVista) del resultado E11: hormigón activo, corte/clip y puntos marcados.

Ejecutar con ▶️ sin argumentos: abre el último estado de ``results/e11_demo/vtk/estado.pvd``
(genéralo antes con scripts/demo_e11_clima_real.py). Requiere PyVista en el .venv:
    .\\.venv\\Scripts\\python.exe -m pip install -r requirements-visor.txt

Argumentos opcionales:
    --tiempo-h 48       estado más cercano a 48 h (por defecto, el último)
    --clip              plano de clip en vez de slice
    --vertida 0         resalta los puntos de todas las vertidas (por defecto: 2)
    --captura DIR       sin ventana: guarda PNG de vista completa, slice y clip en DIR

Cómo mover el plano de corte (modo interactivo, en la ventana):
    * Arrastra con el ratón la FLECHA/NORMAL del plano para desplazarlo.
    * Arrastra el BORDE o el contorno del plano para inclinarlo o girarlo.
    * Arrastra el CÍRCULO del origen para trasladarlo en cualquier dirección.
    * Rueda del ratón = zoom; clic izquierdo y arrastrar sobre el fondo = rotar la vista.

Qué muestra: solo hormigón ACTIVO (el archivo ya no contiene lo no colocado; además se filtra
``mascara_activa == 1``), color = ``temperatura_C`` [°C]. Esferas: núcleo (rojo, grande), voxel próximo
a la cara (verde, T del voxel) y cara exterior con T RECONSTRUIDA (azul, pequeña). La reconstrucción no
es una medición ni la temperatura del voxel. En ParaView: abrir estado.pvd y marcas.pvd; Filters > Slice.
"""
from __future__ import annotations

import argparse
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
PVD_ESTADOS = RAIZ / "results" / "e11_demo" / "vtk" / "estado.pvd"
PVD_MARCAS = RAIZ / "results" / "e11_demo" / "vtk" / "marcas.pvd"
COLORES = {1: "red", 2: "limegreen", 3: "dodgerblue"}
ETIQUETAS = {1: "núcleo", 2: "voxel próximo a la cara", 3: "cara exterior (T reconstruida)"}
RADIOS = {1: 0.28, 2: 0.17, 3: 0.11}          # m: núcleo > voxel próximo > cara (se distinguen aunque estén cerca)
DESPLAZAMIENTO_Z = {1: 0.7, 2: 0.7, 3: -0.7}  # m: separa las etiquetas para que no se tapen


def cargar(pvd: Path, tiempo_h: float | None):
    import pyvista as pv
    lector = pv.PVDReader(str(pvd))
    tiempos = list(lector.time_values)
    elegido = tiempos[-1] if tiempo_h is None else min(tiempos, key=lambda t: abs(t - tiempo_h))
    lector.set_active_time_value(elegido)
    bloque = lector.read()
    malla = bloque[0] if isinstance(bloque, pv.MultiBlock) else bloque
    return malla, float(elegido), tiempos


def solo_activo(malla):
    return malla.threshold(0.5, scalars="mascara_activa") if "mascara_activa" in malla.array_names else malla


def agregar_marcas(plotter, marcas, vertida: int | None):
    import numpy as np
    import pyvista as pv
    etiquetas, puntos = [], []
    for codigo in (1, 2, 3):
        mascara = (marcas["codigo_marca"] == codigo)
        if vertida is not None:
            mascara &= (marcas["vertida_id"] == vertida)
        sub = marcas.extract_points(mascara, adjacent_cells=False)
        if sub.n_points == 0:
            continue
        esferas = sub.glyph(geom=pv.Sphere(radius=RADIOS[codigo]), scale=False, orient=False)
        plotter.add_mesh(esferas, color=COLORES[codigo], label=ETIQUETAS[codigo])
        for xyz, temp in zip(sub.points, sub["temperatura_C"]):
            puntos.append(np.asarray(xyz) + np.array([0.0, 0.0, DESPLAZAMIENTO_Z[codigo]]))
            etiquetas.append(f"{ETIQUETAS[codigo]}: {temp:.1f} °C")
    if puntos:
        plotter.add_point_labels(puntos, etiquetas, font_size=11, shape_opacity=0.65,
                                 always_visible=True, show_points=False)
    plotter.add_legend(bcolor=None, size=(0.45, 0.16), loc="upper left")


def _opciones(malla):
    rango = (float(malla["temperatura_C"].min()), float(malla["temperatura_C"].max()))
    return dict(scalars="temperatura_C", clim=rango, cmap="coolwarm", scalar_bar_args={"title": "T [°C]"})


def crear_visor(malla, marcas, *, vertida=None, clip=False, titulo="", off_screen=False):
    """Plotter con el hormigón activo, el plano de corte interactivo y las marcas."""
    import pyvista as pv
    opciones = _opciones(malla)
    p = pv.Plotter(off_screen=off_screen, window_size=(1200, 850))
    if clip:
        p.add_mesh_clip_plane(malla, normal="y", **opciones)
    else:
        p.add_mesh(malla, opacity=0.12, color="lightgray")
        p.add_mesh_slice(malla, normal="y", **opciones)
    if marcas is not None:
        agregar_marcas(p, marcas, vertida)
    p.add_text(titulo, font_size=10)
    p.add_axes()
    p.camera_position = "iso"
    return p


def capturas(malla, marcas, vertida, titulo, carpeta: Path):
    import pyvista as pv
    carpeta.mkdir(parents=True, exist_ok=True)
    opciones = _opciones(malla)
    centro = malla.center
    vistas = (
        ("3d_completo", lambda p: p.add_mesh(malla, show_edges=False, **opciones)),
        ("3d_slice", lambda p: (p.add_mesh(malla, opacity=0.12, color="lightgray"),
                                p.add_mesh(malla.slice(normal="y", origin=centro), **opciones))),
        ("3d_clip", lambda p: p.add_mesh(
            malla.clip(normal="y", origin=(centro[0], centro[1] + 0.6, centro[2])), **opciones)),
    )
    for nombre, preparar in vistas:
        p = pv.Plotter(off_screen=True, window_size=(1100, 800))
        preparar(p)
        if marcas is not None:
            agregar_marcas(p, marcas, vertida)
        p.add_text(titulo, font_size=10)
        p.add_axes()
        p.camera_position = "iso"
        ruta = carpeta / f"{nombre}.png"
        p.screenshot(str(ruta))
        p.close()
        print(f"Captura: {ruta}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--estados", type=Path, default=PVD_ESTADOS)
    parser.add_argument("--marcas", type=Path, default=PVD_MARCAS)
    parser.add_argument("--tiempo-h", type=float, default=None)
    parser.add_argument("--vertida", type=int, default=2, help="vertida cuyos puntos se resaltan (0 = todas)")
    parser.add_argument("--clip", action="store_true", help="plano de clip interactivo en vez de slice")
    parser.add_argument("--captura", type=Path, default=None, help="carpeta de salida: modo sin ventana")
    args = parser.parse_args()
    try:
        import pyvista  # noqa: F401
    except ImportError:
        print("Falta PyVista. Instálalo con:\n  .\\.venv\\Scripts\\python.exe -m pip install -r requirements-visor.txt")
        return 2
    if not args.estados.is_file():
        print(f"No existe {args.estados}. Ejecuta antes scripts/demo_e11_clima_real.py con ▶️.")
        return 2
    malla, tiempo_h, tiempos = cargar(args.estados, args.tiempo_h)
    malla = solo_activo(malla)
    vertida = args.vertida or None
    marcas = cargar(args.marcas, tiempo_h)[0] if args.marcas.is_file() else None
    titulo = f"t = {tiempo_h:.1f} h · hormigón activo · T [°C] (estados: {len(tiempos)})"
    if args.captura:
        capturas(malla, marcas, vertida, titulo, args.captura)
        return 0
    print(f"Visor abierto: {args.estados} · t = {tiempo_h:.1f} h. Arrastra el plano para cortar; "
          "cierra la ventana para terminar.")
    p = crear_visor(malla, marcas, vertida=vertida, clip=args.clip,
                    titulo=titulo + " · mueve el plano para cortar")
    p.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
