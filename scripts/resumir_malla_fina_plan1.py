"""Compara Plan 1 y PUR de 120 mm en los 24 escenarios refinados."""
import csv
from fractions import Fraction
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATOS = ROOT / "results/cierre_plan1"
SALIDA = ROOT / "reports/Defensa/cierre_plan1/malla_fina_24.csv"
CASOS = (
    ("plan1_336_fino", "60 mm a 336 h", 336, "exacto"),
    ("pur_120mm_mantenido_1008_fino", "120 mm a 1008 h", 1008, "pur_120mm_mantenido_1008"),
)


def resumir(caso, descripcion, horizonte, nominal):
    esperados = {p.name for p in (DATOS / nominal).glob("*.json.gz")}
    archivos = {p.name: p for p in (DATOS / caso).glob("*.json.gz")}
    if len(esperados) != 24 or esperados != archivos.keys():
        raise ValueError(f"{caso}: matriz incompleta ({len(archivos)}/24)")
    filas = []
    for p in sorted(archivos.values()):
        with gzip.open(p, "rt", encoding="utf-8") as f:
            fila = json.load(f)["fila"]
        if (fila["dx_m"] != 0.125 or fila["dt_limite_s"] != 600
                or fila["horizonte_h"] != horizonte or fila["dt_sobre_estable_max"] >= 1):
            raise ValueError(f"Configuración o estabilidad inesperada: {p}")
        filas.append(fila)
    pesos = [Fraction(f["peso"]).limit_denominator(248) for f in filas]
    if sum(pesos) != 1:
        raise ValueError(f"Pesos inválidos: {caso}")
    falla = [f["M_T_c"] > 70 or f["M_delta_c"] > 20 for f in filas]
    riesgo = sum((w for w, excede in zip(pesos, falla) if excede), Fraction())
    mt = max(filas, key=lambda f: f["M_T_c"])
    md = max(filas, key=lambda f: f["M_delta_c"])
    return dict(caso=descripcion, horizonte_h=horizonte, escenarios=24,
                M_T_c=mt["M_T_c"], escenario_T=mt["escenario"], t_T_h=mt["t_M_T_h"],
                M_delta_c=md["M_delta_c"], escenario_delta=md["escenario"], t_delta_h=md["t_M_delta_h"],
                excedencias=sum(falla), probabilidad=str(riesgo),
                riesgo_pct=float(riesgo * 100), cumple=riesgo <= Fraction(1, 20),
                razon_estabilidad_max=max(f["dt_sobre_estable_max"] for f in filas),
                balance_relativo_max=max(f["balance_relativo"] for f in filas))


if __name__ == "__main__":
    filas = [resumir(*caso) for caso in CASOS]
    with SALIDA.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(filas[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(filas)
