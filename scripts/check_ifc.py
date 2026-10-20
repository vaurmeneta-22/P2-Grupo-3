"""Inspecciona las vertidas, propiedades y dimensiones de un archivo IFC.

Al ejecutarse sin argumentos (por ejemplo, con el botón Play de VS Code),
revisa models/ifc/Fundacion_1_con_propiedades.ifc. También acepta una ruta:
    python .\\scripts\\check_ifc.py .\\models\\ifc\\otro_modelo.ifc

Requiere IfcOpenShell y NumPy en el entorno Python activo.
"""
import sys
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np


def main() -> int:
    if len(sys.argv) > 2:
        print(
            "Uso: python .\\scripts\\check_ifc.py [ruta_al_archivo.ifc]",
            file=sys.stderr,
        )
        return 2

    if len(sys.argv) == 2:
        ruta_ifc = Path(sys.argv[1]).expanduser()
    else:
        raiz_proyecto = Path(__file__).resolve().parents[1]
        ruta_ifc = raiz_proyecto / "models" / "ifc" / "Fundacion_1_con_propiedades.ifc"
        print(f"Sin argumento: usaré el IFC predeterminado {ruta_ifc}")

    if not ruta_ifc.is_file():
        print(f"ERROR: no existe el archivo: {ruta_ifc}", file=sys.stderr)
        return 2

    modelo = ifcopenshell.open(str(ruta_ifc))
    vertidas = [
        producto
        for producto in modelo.by_type("IfcProduct")
        if (getattr(producto, "Name", "") or "").startswith("Vertida_")
    ]
    vertidas.sort(key=lambda producto: producto.Name)

    print(f"Archivo: {ruta_ifc}")
    print(f"Esquema: {modelo.schema}")
    print(f"Vertidas encontradas: {len(vertidas)}")
    if not vertidas:
        print("ERROR: no encontré productos con Name='Vertida_XX'.", file=sys.stderr)
        return 1

    configuracion = ifcopenshell.geom.settings()
    configuracion.set("use-world-coords", True)
    minimos_globales = []
    maximos_globales = []
    errores_geometria = 0

    for vertida in vertidas:
        print("\n" + "=" * 68)
        print(f"Nombre IFC: {vertida.Name}")
        print(f"Clase IFC: {vertida.is_a()}")
        print(f"GUID: {vertida.GlobalId}")

        conjuntos = ifcopenshell.util.element.get_psets(vertida)
        if conjuntos:
            print("Conjuntos de propiedades:")
            for nombre, propiedades in conjuntos.items():
                print(f"  {nombre}: {propiedades}")
        else:
            print("Conjuntos de propiedades: ninguno")

        try:
            forma = ifcopenshell.geom.create_shape(configuracion, vertida)
            vertices = np.asarray(forma.geometry.verts, dtype=float).reshape(-1, 3)
            if vertices.size == 0:
                raise ValueError("la geometría no contiene vértices")
            minimo = vertices.min(axis=0)
            maximo = vertices.max(axis=0)
            minimos_globales.append(minimo)
            maximos_globales.append(maximo)
            print(f"Caja mínima [m]: {np.round(minimo, 4)}")
            print(f"Caja máxima [m]: {np.round(maximo, 4)}")
            print(f"Dimensiones [m]: {np.round(maximo - minimo, 4)}")
        except Exception as error:  # informa el elemento y continúa con los demás
            errores_geometria += 1
            print(f"ERROR al leer geometría: {error}")

    if minimos_globales:
        minimo_total = np.min(np.vstack(minimos_globales), axis=0)
        maximo_total = np.max(np.vstack(maximos_globales), axis=0)
        print("\n" + "=" * 68)
        print("RESUMEN GLOBAL")
        print(f"Caja mínima [m]: {np.round(minimo_total, 4)}")
        print(f"Caja máxima [m]: {np.round(maximo_total, 4)}")
        print(f"Dimensiones globales [m]: {np.round(maximo_total - minimo_total, 4)}")

    print(f"\nNúmero de vertidas: {len(vertidas)}")
    print(f"Errores de geometría: {errores_geometria}")
    return 1 if errores_geometria else 0


if __name__ == "__main__":
    raise SystemExit(main())
