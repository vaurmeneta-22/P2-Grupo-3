r"""Inspeccion de SOLO LECTURA del IFC base (no lo modifica ni lo voxeliza).

    python .\scripts\inspeccionar_ifc_e04.py [ruta.ifc]

Informa unidades, vertidas (clase, nombre, GlobalId), propiedades y, por
geometria triangulada: caja envolvente y volumen de malla. La comparacion
volumen de malla / volumen de caja muestra cuanto se aparta cada vertida de
un bloque rectangular (octagonos, troncos y cilindro teselado: <1). NO es una
validacion geometrica completa (no se revisa solape, vacios ni hermeticidad).
"""
from pathlib import Path
import sys

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.unit
import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
IFC_BASE = RAIZ / "models" / "ifc" / "Fundacion_1_con_propiedades.ifc"


def volumen_malla(vertices, caras):
    a, b, c = (vertices[caras[:, i]] for i in range(3))
    return abs(float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)


def main() -> int:
    ruta = Path(sys.argv[1]) if len(sys.argv) > 1 else IFC_BASE
    modelo = ifcopenshell.open(str(ruta))
    escala = ifcopenshell.util.unit.calculate_unit_scale(modelo)
    print(f"Archivo: {ruta}\nEsquema: {modelo.schema}  |  escala de longitud a metros: {escala}")
    for u in modelo.by_type("IfcSiUnit"):
        print(f"  Unidad {u.UnitType}: {u.Prefix or ''}{u.Name}")
    vertidas = sorted((p for p in modelo.by_type("IfcProduct")
                       if (p.Name or "").startswith("Vertida_")), key=lambda p: p.Name)
    clases = sorted({p.is_a() for p in vertidas})
    print(f"Vertidas: {len(vertidas)}  clases: {clases}")
    cfg = ifcopenshell.geom.settings()
    cfg.set("use-world-coords", True)
    problemas = 0
    for p in vertidas:
        forma = ifcopenshell.geom.create_shape(cfg, p).geometry
        v = np.asarray(forma.verts, float).reshape(-1, 3)
        f = np.asarray(forma.faces, int).reshape(-1, 3)
        dim = v.max(0) - v.min(0)
        vol, caja = volumen_malla(v, f), float(np.prod(dim))
        vu = np.unique(np.round(v, 4), axis=0)
        z_niveles = np.unique(vu[:, 2])
        semi = [round(float(np.abs(vu[vu[:, 2] == z][:, :2]).max()), 3) for z in z_niveles]
        problemas += not (np.isfinite(vol) and vol > 0)
        ps = ifcopenshell.util.element.get_psets(p).get("UANDES_HormigonMasivo", {})
        print(f"{p.Name} | {p.is_a()} | {p.GlobalId} | z=[{v[:,2].min():.2f},{v[:,2].max():.2f}] m"
              f" | dim={np.round(dim,3)} m | vol={vol:.3f} m3 | vol/caja={vol/caja:.3f}")
        print(f"    malla: {len(vu)} vertices unicos; semi-ancho por nivel z {[float(z) for z in z_niveles]} = {semi} m")
        print(f"    props: { {k: x for k, x in ps.items() if k != 'id'} }")
    print("\nNOTA: 'geometria procesable' != validacion completa. El solver 2D/3D NO consume este IFC.")
    return 1 if problemas else 0


if __name__ == "__main__":
    raise SystemExit(main())
