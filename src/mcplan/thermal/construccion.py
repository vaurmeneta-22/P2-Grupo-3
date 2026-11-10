"""Conjuntos acumulados de vertidas para las figuras de construcción progresiva.

Los eventos de activación vienen de ``analisis.eventos_activacion(resultado)``
(implementación del equipo); este módulo solo los acumula y los nombra.
"""
from __future__ import annotations


def conjuntos_acumulados(eventos) -> list[tuple[int, ...]]:
    """Conjunto de vertidas colocadas tras cada evento, en orden de colocación."""
    acumulado: list[int] = []
    resultado = []
    for _, ids in eventos:
        acumulado = acumulado + list(ids)
        resultado.append(tuple(acumulado))
    return resultado


def nombre_conjunto(ids) -> str:
    """Texto como ``V1+V2+V3`` en el orden recibido."""
    return "+".join(f"V{i}" for i in ids)
