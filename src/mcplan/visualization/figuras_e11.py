"""Figuras de la evolución térmica E11: generación, temperaturas, flujos y diferencias.

Todas usan el tiempo GLOBAL desde el inicio de obra [h] y salen de una única simulación (observadores
``SeguimientoPuntos``, ``SeguimientoExtremos`` y ``RegistroEnergetico``). Convenciones comunes:

* Color = vertida (V1..V5) en todas las figuras. Aire = gris.
* Temperatura de CELDA (centro del voxel): línea continua (interior) y punteada (voxel próximo a la cara).
* Temperatura RECONSTRUIDA en la cara exterior: línea discontinua. No es la del voxel ni una medición.
* Eventos: línea punteada gris = colocación (+V#); línea rojo discontinua = retiro de protección.
* Flujo: positivo = entra al hormigón, negativo = sale. W/m² es intensidad (por unidad de área);
  W es la potencia integrada sobre las celdas seleccionadas.
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np

COLORES = {1: "#2c7fb8", 2: "#2aa7b5", 3: "#4daf4a", 4: "#e8a33d", 5: "#d7301f"}
COLOR_AIRE = "#7f7f7f"
COLOR_RETIRO = "#c0392b"
ESTILO_INTERIOR, ESTILO_VOXEL, ESTILO_CARA = "-", ":", "--"
ETIQUETA_TIEMPO = "Tiempo desde inicio de obra [h]"


def _color(i: int) -> str:
    return COLORES.get(i, "#555555")


def marcar_eventos(eje, activaciones, retiros, *, con_texto=True, solo_vertida=None):
    """Colocaciones (gris punteado) y retiros de protección (rojo discontinuo) con su hora."""
    for hora, ids in activaciones:
        if solo_vertida is not None and solo_vertida not in ids:
            continue
        eje.axvline(hora, color="#777777", linestyle=":", linewidth=0.9, zorder=1)
        if con_texto:
            eje.annotate("+" + "+".join(f"V{i}" for i in ids), xy=(hora, 1.0), xycoords=("data", "axes fraction"),
                         xytext=(2, -3), textcoords="offset points", rotation=90, va="top", ha="left",
                         fontsize=7, color="#555555")
    for hora, nombres in retiros:
        if solo_vertida is not None and not any(n.startswith(f"V{solo_vertida}:") for n in nombres):
            continue
        eje.axvline(hora, color=COLOR_RETIRO, linestyle=(0, (5, 2)), linewidth=1.0, zorder=1)
        if con_texto:
            eje.annotate("retiro " + ", ".join(nombres), xy=(hora, 0.0), xycoords=("data", "axes fraction"),
                         xytext=(2, 3), textcoords="offset points", rotation=90, va="bottom", ha="left",
                         fontsize=7, color=COLOR_RETIRO)


def _ajustar(eje, horas_final, ylabel=None, titulo=None):
    eje.set_xlim(0.0, horas_final)
    if ylabel:
        eje.set_ylabel(ylabel)
    if titulo:
        eje.set_title(titulo, fontsize=10)
    eje.grid(True, alpha=0.3)


def figura_hidratacion(reg, datos, horas_final, ruta):
    """Qdot_h(t) por vertida [W/m³] y potencia total generada [kW]."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True)
    horas = reg.horas
    for i in reg.ids:
        a1.plot(horas, reg.qdot_w_m3[i], color=_color(i), linewidth=1.7, label=f"V{i}")
        pico = reg.pico(reg.qdot_w_m3[i])
        if pico:
            a1.plot(*pico, marker="o", color=_color(i), markeredgecolor="black", zorder=3)
    marcar_eventos(a1, datos["activaciones"], [], con_texto=True)
    _ajustar(a1, horas_final, "Qdot_h [W/m³]", "Generación volumétrica por vertida (pulso de hidratación según la EDAD de cada vertida)")
    a1.legend(ncol=5, fontsize=8, loc="upper right")
    apilado = [np.nan_to_num(np.asarray(reg.potencia_vertida_w[i]) / 1000.0) for i in reg.ids]
    a2.stackplot(horas, apilado, colors=[_color(i) for i in reg.ids], alpha=0.55,
                 labels=[f"V{i}" for i in reg.ids])
    a2.plot(horas, reg.potencia_total_kw(), color="black", linewidth=1.5, label="Total del conjunto activo")
    pico = reg.pico(reg.potencia_total_kw())
    if pico:
        a2.plot(*pico, marker="o", color="black", markerfacecolor="yellow", zorder=4)
        a2.annotate(f"máx {pico[1]:.0f} kW a las {pico[0]:.1f} h", xy=pico, xytext=(10, -22),
                    textcoords="offset points", fontsize=8, arrowprops={"arrowstyle": "-"})
    marcar_eventos(a2, datos["activaciones"], [], con_texto=False)
    _ajustar(a2, horas_final, "Potencia generada [kW]", "Potencia total generada: Σ Qdot_h · V_voxel (apilada por vertida)")
    a2.set_xlabel(ETIQUETA_TIEMPO)
    a2.legend(ncol=6, fontsize=8, loc="upper right")
    fig.suptitle("Generación de calor por hidratación", fontsize=13)
    fig.text(0.01, 0.005, "Cada vertida repite el mismo pulso desplazado a su hora de colocación (líneas punteadas). "
             "Es generación de CALOR (potencia), no de temperatura.", fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(ruta, dpi=140)
    return fig


def figura_temperaturas_por_vertida(seguimiento, datos, horas_final, ruta):
    """Un panel por vertida: interior, voxel próximo, cara exterior reconstruida y aire."""
    puntos = seguimiento.puntos
    fig, ejes = plt.subplots(3, 2, figsize=(13, 11), sharex=True, sharey=True)
    plano = ejes.ravel()
    t_aire = np.asarray(datos["aire_horas"])
    for eje, punto in zip(plano, puntos):
        s = seguimiento.serie(punto.vertida_id)
        h = s["tiempos_s"] / 3600.0
        eje.plot(t_aire, datos["aire_c"], color=COLOR_AIRE, linewidth=1.0, label="Aire (snapshot real)")
        c = _color(punto.vertida_id)
        eje.plot(h, s["interior_c"], color=c, linestyle=ESTILO_INTERIOR, linewidth=1.8, label="Interior (celda)")
        eje.plot(h, s["celda_proxima_c"], color=c, linestyle=ESTILO_VOXEL, linewidth=1.3,
                 label="Voxel próximo a la cara (celda)")
        eje.plot(h, s["cara_estimada_c"], color="black", linestyle=ESTILO_CARA, linewidth=1.2,
                 label=f"Cara {punto.cara} RECONSTRUIDA")
        marcar_eventos(eje, datos["activaciones"], datos["retiros"], con_texto=False, solo_vertida=punto.vertida_id)
        _ajustar(eje, horas_final, "T [°C]", f"Vertida {punto.vertida_id}")
    for eje in plano[len(puntos):]:
        eje.axis("off")
    manijas, etiquetas = plano[0].get_legend_handles_labels()
    plano[-1].legend(manijas, etiquetas, loc="center", fontsize=9, title="Curvas (todas en °C)")
    for eje in ejes[-1]:
        eje.set_xlabel(ETIQUETA_TIEMPO)
    fig.suptitle("Temperaturas por vertida: celdas vs. cara reconstruida vs. aire real "
                 "(sin curva antes de colocarse; punteada gris = colocación, roja = retiro)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ruta, dpi=140)
    return fig


def figura_conjunto(extremos, datos, horas_final, ruta):
    """Tmin, Tmedia y Tmax de las celdas activas del conjunto construido."""
    s = extremos.serie("global")
    h = s["tiempos_s"] / 3600.0
    fig, eje = plt.subplots(figsize=(12, 6.2))
    eje.fill_between(h, s["tmin"], s["tmax"], color="#9ecae1", alpha=0.35, linewidth=0,
                     label="Banda Tmin–Tmax de CELDAS activas")
    eje.plot(h, s["tmax"], color="#d62728", linewidth=1.5, label="Tmax (celdas)")
    eje.plot(h, s["tmedia"], color="black", linewidth=1.8, label="Tmedia (celdas)")
    eje.plot(h, s["tmin"], color="#1f77b4", linewidth=1.5, label="Tmin (celdas)")
    eje.plot(datos["aire_horas"], datos["aire_c"], color=COLOR_AIRE, linewidth=1.0, label="Aire (snapshot real)")
    marcar_eventos(eje, datos["activaciones"], datos["retiros"])
    if np.any(np.isfinite(s["tmax"])):
        n = int(np.nanargmax(s["tmax"]))
        eje.plot(h[n], s["tmax"][n], marker="o", color="#d62728", markeredgecolor="black", zorder=4)
        eje.annotate(f"Tmax {s['tmax'][n]:.2f} °C a las {h[n]:.1f} h", xy=(h[n], s["tmax"][n]),
                     xytext=(12, 12), textcoords="offset points", fontsize=8, arrowprops={"arrowstyle": "-"})
    _ajustar(eje, horas_final, "T [°C]", "Conjunto construido: solo hormigón ya colocado")
    eje.set_xlabel(ETIQUETA_TIEMPO)
    eje.legend(fontsize=8, ncol=3, loc="lower right")
    fig.text(0.01, 0.005, "Los saltos de Tmin/Tmedia al entrar una vertida cambian el conjunto considerado; no son error del solver.",
             fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(ruta, dpi=140)
    return fig


def _etiqueta_cara(reg, v, cara, protecciones):
    nombres = [p.nombre for p in protecciones if p.vertida_id == v and cara in p.caras]
    estado = "protegida: " + "+".join(nombres) if nombres else "expuesta"
    return f"V{v} {cara} ({estado})"


def figura_flujos(reg, datos, horas_final, ruta):
    """Flujo por caras representativas (W/m² y W), suelo en serie separada y línea de tiempo de protecciones."""
    fig, ejes = plt.subplots(5, 1, figsize=(12, 17), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2.4, 1.8, 1.8, 1.5]})
    a_q, a_w, a_sq, a_sw, a_p = ejes
    h = reg.horas
    for clave in reg.caras:
        v, cara = clave
        etiqueta = _etiqueta_cara(reg, v, cara, datos["protecciones"])
        a_q.plot(h, reg.cara_w_m2[clave], color=_color(v), linewidth=1.5, label=etiqueta)
        a_w.plot(h, reg.cara_w[clave], color=_color(v), linewidth=1.5, label=etiqueta)
    for eje in (a_q, a_w, a_sq, a_sw, a_p):
        marcar_eventos(eje, datos["activaciones"], datos["retiros"], con_texto=(eje is a_q))
    for eje in (a_q, a_sq):
        eje.axhline(0.0, color="#555555", linewidth=0.8)
    _ajustar(a_q, horas_final, "Flujo por área q [W/m²]\n(+ entra, − sale)",
             "Flujo por unidad de área en la cara lateral seleccionada de cada vertida (promedio de sus celdas expuestas)")
    a_q.legend(fontsize=7, loc="lower right", ncol=1)
    a_w.axhline(0.0, color="#555555", linewidth=0.8)
    _ajustar(a_w, horas_final, "Potencia integrada [W]\n(+ entra, − sale)",
             "Potencia integrada sobre esas mismas celdas (depende de cuántas celdas expuestas tiene cada vertida)")
    if reg.celdas_suelo is not None:
        a_sq.plot(h, reg.suelo_w_m2, color="#8c564b", linewidth=1.6, label="Suelo (base física, cara z−)")
        a_sw.plot(h, reg.suelo_w, color="#8c564b", linewidth=1.6, label="Suelo (base física, cara z−)")
        pico = reg.pico(reg.suelo_w_m2, signo=-1.0)
        if pico:
            a_sq.plot(*pico, marker="o", color="#8c564b", markeredgecolor="black", zorder=3)
            a_sq.annotate(f"máx salida {pico[1]:.1f} W/m² a las {pico[0]:.1f} h", xy=pico, xytext=(10, 10),
                          textcoords="offset points", fontsize=8, arrowprops={"arrowstyle": "-"})
    _ajustar(a_sq, horas_final, "Suelo q [W/m²]\n(+ entra, − sale)", "Intercambio con el suelo (serie separada): flujo por área")
    _ajustar(a_sw, horas_final, "Suelo [W]\n(+ entra, − sale)", "Intercambio con el suelo (serie separada): potencia integrada en la base")
    # línea de tiempo de protecciones
    filas = []
    for p in datos["protecciones"]:
        inicio = datos["inicio_h"][p.vertida_id]
        fin = horas_final if p.retiro_s is None else p.retiro_s / 3600.0
        retiro_h = None if p.retiro_s is None else p.retiro_s / 3600.0
        filas.append((f"V{p.vertida_id} {p.nombre} (caras {', '.join(p.caras)})", inicio, min(fin, horas_final),
                      p.vertida_id, retiro_h))
    for n, (nombre, ini, fin, v, retiro_h) in enumerate(filas):
        a_p.barh(n, fin - ini, left=ini, color=_color(v), alpha=0.65, height=0.6)
        if retiro_h is None:
            nota = " · no se retira"
        elif retiro_h <= horas_final:
            nota = f" · se retira a las {retiro_h:.0f} h"
        else:
            nota = f" · se retiraría a las {retiro_h:.0f} h (fuera de la ventana)"
        a_p.text(ini + 0.5, n, nombre + nota, va="center", fontsize=7)
    a_p.set_yticks([])
    a_p.invert_yaxis()
    _ajustar(a_p, horas_final, None, "Protecciones vigentes (barra = puesta desde la colocación hasta su retiro)")
    a_p.set_xlabel(ETIQUETA_TIEMPO)
    fig.suptitle("Flujos hacia el entorno (convenio: positivo = entra al hormigón, negativo = sale)", fontsize=12)
    fig.text(0.01, 0.003, "Las vertidas se colocan a horas distintas y con temperatura de colocación distinta: a igual hora global "
             "tienen edades y temperaturas distintas (ver flujos_por_edad.png).", fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.012, 1, 0.975))
    fig.savefig(ruta, dpi=130)
    return fig


def figura_por_edad(reg, datos, ruta):
    """Qdot y flujo por área en función de la EDAD de cada vertida (horas desde su colocación)."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8.4), sharex=True)
    h = reg.horas
    for i in reg.ids:
        edad = h - datos["inicio_h"][i]
        sel = edad >= 0
        a1.plot(edad[sel], np.asarray(reg.qdot_w_m3[i])[sel], color=_color(i), linewidth=1.6, label=f"V{i}")
    for clave in reg.caras:
        v, cara = clave
        edad = h - datos["inicio_h"][v]
        sel = (edad >= 0) & np.isfinite(reg.cara_w_m2[clave])
        a2.plot(edad[sel], np.asarray(reg.cara_w_m2[clave])[sel], color=_color(v), linewidth=1.5,
                label=_etiqueta_cara(reg, v, cara, datos["protecciones"]))
    a1.set_ylabel("Qdot_h [W/m³]")
    a1.set_title("Generación por hidratación vs. edad: coinciden entre vertidas (mismo modelo de pulso)", fontsize=10)
    a2.axhline(0.0, color="#555555", linewidth=0.8)
    a2.set_ylabel("q cara lateral [W/m²]\n(+ entra, − sale)")
    a2.set_title("Flujo por área vs. edad: difiere por protección, temperatura de colocación y aire a distinta hora", fontsize=10)
    a2.set_xlabel("Edad de cada vertida [h desde su colocación]")
    for eje in (a1, a2):
        eje.grid(True, alpha=0.3)
    a1.legend(ncol=5, fontsize=8)
    a2.legend(fontsize=7, loc="lower right")
    fig.suptitle("Comparación a igual edad (no a igual hora global)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ruta, dpi=140)
    return fig


def figura_diferencias(seguimiento, extremos, datos, horas_final, ruta):
    """ΔT interior−superficie por vertida y diferencias del dominio, con sus máximos y horas."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11.5, 9), sharex=True)
    for k, fila in enumerate(seguimiento.resumen()):
        p = fila["punto"]
        s = seguimiento.serie(p.vertida_id)
        h = s["tiempos_s"] / 3600.0
        a1.plot(h, s["dt_cara_c"], color=_color(p.vertida_id), linewidth=1.5, label=f"V{p.vertida_id}")
        if fila["dt_cara"] is not None:
            r = fila["dt_cara"]
            a1.plot(r["tiempo_h"], r["diferencia_con_signo_c"], marker="o", color=_color(p.vertida_id),
                    markeredgecolor="black", zorder=3)
            # etiquetas escalonadas: varias vertidas pueden tener su máximo a la misma hora
            a1.annotate(f"V{p.vertida_id}: máx {r['max_abs_c']:.1f} °C a las {r['tiempo_h']:.0f} h",
                        xy=(r["tiempo_h"], r["diferencia_con_signo_c"]), xytext=(-165, 14 + 13 * k),
                        textcoords="offset points", fontsize=7, color=_color(p.vertida_id),
                        arrowprops={"arrowstyle": "-", "color": _color(p.vertida_id), "linewidth": 0.6})
    marcar_eventos(a1, datos["activaciones"], datos["retiros"], con_texto=False)
    a1.axhline(0.0, color="#555555", linewidth=0.8)
    _ajustar(a1, horas_final, "ΔT = T interior − T cara\nreconstruida [°C]",
             "Diferencia interior − superficie por vertida (círculo = máx |ΔT|, con su hora)")
    a1.legend(ncol=5, fontsize=8)
    s = extremos.serie("global")
    h = s["tiempos_s"] / 3600.0
    for campo, etiqueta, color in (("dt_global_c", "Tmax − Tmin entre CELDAS activas", "#d62728"),
                                   ("contraste_abs_c", "Máx |T celda − T cara exterior reconstruida|", "#1f77b4")):
        y = s[campo]
        a2.plot(h, y, color=color, linewidth=1.6, label=etiqueta)
        if np.any(np.isfinite(y)):
            n = int(np.nanargmax(y))
            a2.plot(h[n], y[n], marker="o", color=color, markeredgecolor="black", zorder=3)
            a2.annotate(f"máx {y[n]:.2f} °C a las {h[n]:.1f} h", xy=(h[n], y[n]), xytext=(12, -14 if campo == "dt_global_c" else 10),
                        textcoords="offset points", fontsize=8, color=color, arrowprops={"arrowstyle": "-", "color": color})
    marcar_eventos(a2, datos["activaciones"], datos["retiros"], con_texto=False)
    _ajustar(a2, horas_final, "Diferencia [°C]", "Dominio activo completo")
    a2.set_xlabel(ETIQUETA_TIEMPO)
    a2.legend(fontsize=8, loc="lower right")
    fig.suptitle("Diferencias térmicas (en °C entre puntos; no son gradientes ni criterio de fisuración)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ruta, dpi=140)
    return fig
