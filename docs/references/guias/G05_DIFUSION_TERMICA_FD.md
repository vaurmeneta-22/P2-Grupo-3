# Guía G05 — Difusión térmica y solver cartesiano

## Ecuación

Para propiedades constantes:

\[
\rho c_p\frac{\partial T}{\partial t}
=
k\nabla^2T+\dot Q.
\]

También:

\[
\frac{\partial T}{\partial t}
=
\alpha\nabla^2T+
\frac{\dot Q}{\rho c_p},
\qquad
\alpha=\frac{k}{\rho c_p}.
\]

## Interpretación por balance

Para una celda:

\[
C_i\frac{dT_i}{dt}
=
\sum_jG_{ij}(T_j-T_i)+\dot Q_iV_i+Q_{\mathrm{borde},i},
\]

con:

\[
C_i=\rho_i c_{p,i}V_i
\]

y, para una interfaz simple,

\[
G_{ij}\approx\frac{k_{ij}A_{ij}}{d_{ij}}.
\]

Esta forma ayuda a evitar errores conceptuales cuando:

- cambian materiales;
- aparece suelo;
- hay superficies;
- se agregan resistencias térmicas.

## Condiciones de borde

### Dirichlet

\[
T=T_b
\]

Temperatura impuesta.

### Neumann

\[
-k\nabla T\cdot n=q_n
\]

Flujo impuesto. Caso adiabático:

\[
q_n=0.
\]

### Robin

\[
-k\nabla T\cdot n=h(T-T_\infty)
\]

Convección.

## Fuente versus condición de borde

### Fuente volumétrica

Calor de hidratación:

\[
\dot Q_h\quad[\mathrm{W/m^3}]
\]

o unidades equivalentes coherentes.

### Flujo superficial

Radiación solar absorbida:

\[
q''_{\mathrm{solar}}\quad[\mathrm{W/m^2}]
\]

### Convección

No es una “fuente interna”: es transferencia a través de la frontera.

## Estabilidad explícita

Para difusión pura, grilla uniforme e isotrópica:

### 1D

\[
Fo=\frac{\alpha\Delta t}{\Delta x^2}\le\frac12.
\]

### 2D, \(\Delta x=\Delta y\)

\[
Fo\le\frac14.
\]

### 3D, \(\Delta x=\Delta y=\Delta z\)

\[
Fo\le\frac16.
\]

Forma más general:

\[
\alpha\Delta t
\left(
\frac1{\Delta x^2}+
\frac1{\Delta y^2}+
\frac1{\Delta z^2}
\right)
\le\frac12.
\]

Para balances con propiedades heterogéneas y conductancias, una guía local útil es:

\[
\Delta t
\lesssim
\frac{C_i}{\sum_jG_{ij}+G_{\mathrm{borde},i}}.
\]

La implementación debe verificar estabilidad empírica y convergencia; no basta con copiar
una fórmula.

## Tests recomendados

1. solución sinusoidal transiente 1D;
2. solución sinusoidal transiente 2D;
3. solución sinusoidal transiente 3D;
4. sistema adiabático sin fuente: conservación de energía;
5. fuente uniforme adiabática: aumento lineal de temperatura;
6. condición de flujo;
7. condición convectiva contra solución 1D conocida.

## Método de Schmidt

El Anexo 4 de la especificación ICH describe un método simplificado de diferencias finitas
para hormigón masivo y señala que puede aplicarse en 1D, 2D y 3D e incorporar secuencia
constructiva, aislación e interacción con otros materiales.

Esa referencia justifica bien la filosofía del proyecto, aunque nuestro solver puede usar
una formulación de balance algo más general.
