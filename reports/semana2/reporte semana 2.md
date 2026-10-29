# Reporte Semana 2 — Solver térmico 2D/3D y verificación

**Proyecto:** Proyecto 2 — Planificación de vertidas de hormigón masivo  
**Equipo:** Grupo 3  
**Integrantes:** Monserrat Cubillos Tome, Josefina Ignacia Muro Teuber y Vicente Andrés Urmeneta Campos  
**Repositorio:** `P2-Grupo-3`  
**Fecha de registro:** 29 de octubre de 2026

## 1. Propósito y alcance

Durante la Semana 2 el proyecto pasó de un ejemplo didáctico unidimensional a un solver explícito de conducción térmica en dominios cartesianos 2D y 3D. Se incorporaron condiciones de borde seleccionables por cara, generación volumétrica de calor, materiales con propiedades distintas y una batería automatizada de verificaciones numéricas.

Los demos usan bloques rectangulares idealizados. El modelo IFC de la fundación se conserva como modelo BIM de referencia, pero todavía no se convierte en una malla térmica ni se acopla al solver. Por ello, los resultados de estos demos verifican el método y el código; no son una predicción térmica calibrada de la fundación real.

## 2. Archivos y responsabilidades

| Ruta | Contenido |
| --- | --- |
| `src/mcplan/thermal/diffusion_nd.py` | Solver explícito por balance de energía en volúmenes/celdas cartesianas 2D y 3D. |
| `src/mcplan/thermal/verificacion.py` | Soluciones de referencia y verificaciones para modos sinusoidales y fuente uniforme. |
| `scripts/demo_difusion_2d.py` | Demo 2D configurable; en la versión de trabajo documentada modela un bloque de 4 × 2 m y una fuente temporal. |
| `scripts/demo_difusion_3d.py` | Demo 3D con casos base y multidimensionales. |
| `scripts/verificar_E05.py` | Ejecuta las verificaciones de solución exacta 2D/3D y fuente uniforme adiabática, e informa errores y tolerancias. |
| `scripts/run_tests.py` | Descubre y ejecuta las pruebas `test_*.py` de `tests/`, pensado también para abrirse y ejecutarse desde VS Code. |
| `scripts/inspeccionar_ifc_e04.py` | Inspecciona en modo de solo lectura unidades, vertidas, propiedades y geometría triangulada del IFC. |
| `tests/test_diffusion_1d.py` | Pruebas del solver 1D, sus fuentes y sus bordes. |
| `tests/test_diffusion_nd.py` | Pruebas del solver 2D/3D, solución de referencia, bordes, fuentes, materiales y estabilidad. |
| `models/ifc/Fundacion_1_con_propiedades.ifc` | Modelo IFC4 de la fundación con cinco vertidas y propiedades del curso. |
| `models/bim/Fundacion_Grupo3_Presentacion.blend` | Escena Blender/Bonsai de presentación asociada al trabajo BIM. |
| `ejercicios_locales/fuente_variable_2d/` | Ejercicio independiente agregado desde GitHub sobre una fuente variable y el cálculo de un `Qmax` límite. |
| `ejercicios_locales/fuente_variable_3d/` | Ejercicio independiente de clase: fuente temporal en una caja 3D y comparación de condiciones de borde. |

## 3. Balance de energía y discretización

El solver multidimensional actualiza la temperatura de cada celda mediante un balance de energía explícito:

\[
\rho_i c_{p,i} V_i\frac{dT_i}{dt}
=\sum_{\text{caras}}\Phi_{i,\text{cara}}+\dot Q_i V_i.
\]

Aquí \(\Phi\) es la potencia que entra a la celda por una cara, en W; \(\dot Q\) es la generación volumétrica, en W/m³; y \(V_i\) es el volumen representado por la celda. En 2D se define además una profundidad para que el balance tenga volumen y energía en unidades físicas.

Para dos celdas vecinas, las resistencias conductivas de media celda se suman. Esto equivale a usar una conductividad efectiva de media armónica en la interfaz. El flujo que sale de una celda entra con el mismo valor y signo opuesto a la celda vecina, por lo que el intercambio interno conserva energía salvo redondeo numérico. Las propiedades térmicas pueden ser escalares o arreglos por celda, lo que permite representar más de un material.

El avance temporal es Euler explícito. El límite estable se calcula a partir de la capacidad térmica de cada celda y la suma de conductancias que la conectan:

\[
\Delta t\leq\min_i\frac{C_i}{\sum_j G_{ij}},\qquad C_i=\rho_i c_{p,i}V_i.
\]

El límite incorpora los bordes Dirichlet y Robin y las interfaces entre materiales. Así, una cara convectiva intensa o una celda con capacidades distintas puede exigir un paso menor que el estimado para un interior homogéneo.

## 4. Condiciones de borde

En 2D se configura cada una de las cuatro caras; en 3D, cada una de las seis. La convención de Neumann multidimensional es flujo superficial entrante \(q_{in}\) en W/m²: un valor positivo calienta el dominio, uno negativo lo enfría y cero representa una cara adiabática.

| Condición | Interpretación en el solver multidimensional |
| --- | --- |
| Dirichlet | Temperatura prescrita en la cara. El intercambio se calcula con la resistencia conductiva entre el centro de la celda y la superficie. |
| Neumann | Flujo superficial entrante impuesto; la potencia en la cara es \(q_{in}A\). |
| Neumann cero / aislado | No cruza calor por esa cara; condición adiabática. |
| Robin | Convección con temperatura ambiente y coeficiente \(h\); conducción en media celda y convección se representan como resistencias en serie. |

La convención Neumann depende del solver: el módulo 1D recibe un gradiente \(dT/dx\) en K/m, mientras que el solver 2D/3D recibe directamente el flujo superficial en W/m². La conversión depende de la orientación de la cara y de la ley de Fourier.

## 5. Fuentes de calor

El solver multidimensional acepta una fuente volumétrica uniforme, distinta por celda o función del tiempo. Al ser volumétrica, su aporte a una celda es \(\dot Q V\) en W. Esto se distingue de un flujo prescrito en una superficie, cuyo aporte es \(q''A\) en W.

En el demo 2D de la versión de trabajo registrada se usa un bloque de 4 × 2 m, con temperatura inicial de 20 °C y Dirichlet de 20 °C en sus cuatro caras. La fuente editable implementa:

\[
\dot Q(t)=
\begin{cases}
\dfrac{Q_{max}}{2}\left(1-\cos\left(\dfrac{\pi t}{T}\right)\right),&0\le t\le 2T,\\
0,&t>2T.
\end{cases}
\]

En esa configuración \(Q_{max}=1300\) W/m³ y \(T=24\) h: la fuente alcanza su máximo al día 1, retorna a cero al día 2 y permanece apagada. El horizonte configurado para este demo es 50 h. Estas cifras son parámetros ilustrativos del caso de estudio, no una tasa de hidratación calibrada.

### Ejercicio independiente agregado desde GitHub

`ejercicios_locales/fuente_variable_2d/` contiene un segundo ejercicio, separado de los demos principales y del módulo `src/mcplan/`. Usa su propia copia del solver y plantea encontrar el mayor \(Q_{max}\) que mantenga el centro bajo 70 °C. Su pulso se expresa con duración \(\tau=48\) h:

\[
\dot Q(t)=\frac{Q_{max}}{2}\left(1-\cos\left(\frac{2\pi t}{\tau}\right)\right),\quad 0\le t\le\tau,
\qquad \dot Q(t)=0\quad\text{fuera de ese intervalo}.
\]

El código estima \(Q_{max}\) por linealidad, lo vuelve a comprobar mediante una corrida y contempla sensibilidad separada a la malla y al paso temporal. Incluye un README, un solver propio, un script y pruebas locales. Sus pruebas locales **no forman parte** de las 46 pruebas de la batería principal descrita en la sección 7 y no se registran aquí como ejecutadas.

### Ejercicio independiente 3D agregado desde GitHub

El commit `4b669a2` incorporó `ejercicios_locales/fuente_variable_3d/`, un ejercicio de clase separado del solver principal y de E04/E05. Modela una caja rectangular de 3 × 3 × 0,60 m, inicialmente a 20 °C, con una fuente uniforme en el espacio y variable en el tiempo. Compara Dirichlet a 20 °C, Dirichlet a 30 °C, Neumann con flujo cero y Robin con distintos valores exploratorios de `beta = h/k`.

La fuente usa un pulso cosenoidal de duración editable `τ = 48 h`; el `Qmax` de 1282,565913 W/m³ se toma del ejercicio 2D y no se recalibra para esta caja 3D. El script genera series temporales, cortes espaciales y animaciones GIF, además de un resumen comparativo y un análisis de sensibilidad. Incluye una copia propia de `solver_ejercicio.py`, por lo que sus cambios y pruebas son independientes de `src/mcplan/thermal/diffusion_nd.py` y de las 46 pruebas principales. Según la confirmación del equipo, las **9 pruebas locales del ejercicio 3D se ejecutaron y todas aprobaron**; no se mezclan con los resultados de la batería principal.

## 6. Verificación con soluciones de referencia

Para una placa 2D con Dirichlet uniforme \(T_0\), fuente nula y condición inicial sinusoidal, la referencia implementada es:

\[
T(x,y,t)=T_0+A\sin\left(\frac{\pi x}{L_x}\right)
\sin\left(\frac{\pi y}{L_y}\right)e^{-\alpha\lambda t},
\quad
\lambda=\left(\frac{\pi}{L_x}\right)^2+\left(\frac{\pi}{L_y}\right)^2.
\]

En 3D se añade \(\sin(\pi z/L_z)\) y el término correspondiente a \(z\) en \(\lambda\). La prueba compara todo el campo numérico con la solución analítica y estima una tolerancia que considera los errores espacial y temporal.

Para un dominio adiabático con fuente uniforme constante, la temperatura de referencia aumenta a razón de \(\dot Q/(\rho c_p)\), y el cambio de energía es \(\Delta E=\dot Q Vt\). Estas referencias comprueban tanto el campo térmico como el balance de energía.

## 7. Resultados de la batería ejecutada

El 29 de octubre de 2026 se ejecutó desde la raíz del repositorio, antes de integrar el commit `4b669a2`:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_tests.py
```

**Resultado:** 46 pruebas ejecutadas, 46 aprobadas (`OK`). La batería cubrió, entre otros casos, solución exacta y convergencia 1D, Dirichlet, Neumann, Robin, estabilidad, fuente, conservación adiabática, soluciones de referencia 2D/3D y transferencia entre materiales.

También se ejecutó el verificador específico E05:

```powershell
.\.venv\Scripts\python.exe .\scripts\verificar_E05.py
```

| Caso verificado | Error absoluto máximo | Tolerancia | Resultado |
| --- | ---: | ---: | --- |
| Solución sinusoidal 2D | 0.0086334 °C | 0.039340 °C | Aprobado |
| Solución sinusoidal 3D | 0.017351 °C | 0.083606 °C | Aprobado |
| Fuente uniforme adiabática | 1.42 × 10⁻¹⁴ °C | 2.50 × 10⁻⁹ °C | Aprobado |

En el caso de fuente uniforme, el cambio de energía numérico fue \(1.296000\times10^6\) J, igual al esperado con \(\dot Q Vt\) a la precisión reportada. El error relativo de energía fue \(4.94\times10^{-16}\), menor que la tolerancia de \(10^{-10}\). El verificador terminó con `TODAS LAS VERIFICACIONES APROBADAS`.

## 8. Correspondencia con los casos técnicos E04/E05

| Caso o capacidad | Código y evidencia disponible |
| --- | --- |
| Solver 2D | `diffusion_nd.py`, `demo_difusion_2d.py` y prueba sinusoidal 2D. |
| Solver 3D | `diffusion_nd.py`, `demo_difusion_3d.py` y prueba sinusoidal 3D. |
| Dirichlet | Configurable por cara; hay pruebas de bordes y solución exacta con Dirichlet. |
| Neumann | Flujo entrante configurable, incluido flujo cero adiabático; pruebas de signo y balance. |
| Robin | Intercambio convectivo con resistencias en serie; pruebas de dirección, estado estacionario y caso \(h=0\). |
| Fuente volumétrica | Fuente uniforme y por celda; verificación E05 de energía y temperatura. |
| Conservación adiabática | Prueba de conservación de energía y comparación con el balance de fuente. |
| Dos materiales | Propiedades por celda, conductancia de interfaz y prueba de flujo continuo/resistencia. |
| Modelo BIM base | IFC4 en `models/ifc/` y escena de presentación en `models/bim/`; el solver no los voxeliza todavía. |
| Ejercicio adicional de fuente variable 3D | `ejercicios_locales/fuente_variable_3d/`; actividad independiente, con solver y pruebas locales propios, incorporada por `4b669a2`. No forma parte de la batería principal. |

## 9. Modelo BIM/IFC de referencia

El repositorio conserva `models/ifc/Fundacion_1_con_propiedades.ifc`, con las cinco vertidas `IfcFooting` y el conjunto de propiedades `UANDES_HormigonMasivo` descritos en el reporte de Semana 1. El modelo visual de Blender/Bonsai se conserva en `models/bim/Fundacion_Grupo3_Presentacion.blend`.

Se agregó `scripts/inspeccionar_ifc_e04.py` como inspección IFC de solo lectura. Informa unidades, elementos identificables, propiedades, caja envolvente y volumen de la malla triangulada. La propia documentación del script aclara que esta inspección no constituye una validación geométrica completa. El IFC y la escena `.blend` son entregables relacionados pero distintos; modificar uno no actualiza automáticamente el otro.

## 10. Cómo repetir las verificaciones

Con el entorno del proyecto activo, desde la raíz:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_tests.py
.\.venv\Scripts\python.exe .\scripts\verificar_E05.py
```

Para inspeccionar el IFC base:

```powershell
.\.venv\Scripts\python.exe .\scripts\inspeccionar_ifc_e04.py
```

Los ejercicios independientes de fuente variable 2D y 3D tienen instrucciones de ejecución y pruebas propias en sus respectivos README; no se incluyen en la corrida principal de 46 pruebas. El equipo reporta que las 9 pruebas locales del ejercicio 3D aprobaron.

## 11. Límites y siguiente integración técnica

- Los demos 2D/3D usan bloques rectangulares cartesianos, no la geometría de la fundación IFC.
- No existe todavía un flujo automatizado IFC → voxelización → campos térmicos por vertida.
- No se incorporaron interfaces BIM con suelo, moldajes o aislación en el dominio térmico principal.
- La fuente temporal de los demos es una forma impuesta para estudiar el solver; no equivale a una ley calibrada de hidratación del hormigón.
- El ejercicio de \(Q_{max}\) vive en una carpeta independiente y mantiene una copia propia del solver; sus cambios no se propagan al módulo principal.

## 12. Estado de Git al registrar las pruebas

La batería principal que se reporta se ejecutó sobre el árbol de trabajo con `HEAD` en `8d0b4b1` (`Agrega ejercicio 2D con fuente variable y cálculo del máximo admisible`), antes de integrar el commit remoto `4b669a2` (`Agrega ejercicio térmico 3D con fuente variable y comparación de condiciones de borde`). La rama local se actualizó por avance rápido hasta `4b669a2`; no se creó un commit nuevo. Se conservaron las modificaciones locales en el IFC, los demos 1D/2D y el solver 2D/3D, además de `scripts/run_tests.py` y este reporte sin seguimiento. El SHA del commit remoto no identifica por sí solo todas estas modificaciones. Antes de declarar el commit evaluado, hay que incorporar al commit los archivos que el equipo decida entregar.

