# Ejercicio 2D con fuente variable en el tiempo

Ejercicio de clases **independiente** de las entregas E04/E05: no modifica ni depende de `src/mcplan/`,
`tests/`, los demos ni `scripts/verificar_E05.py`. Usa su propia copia del solver (`solver_ejercicio.py`).

## Problema físico

Rectángulo de **4 m (x) × 2 m (y)**, con profundidad de 1 m fuera del plano (las energías salen en J por
esa profundidad). Temperatura inicial uniforme de **20 °C** y los cuatro bordes con **Dirichlet a 20 °C**.
Material homogéneo de propiedades constantes (los mismos valores del demo 2D del repositorio):

| Propiedad | Valor | Unidad |
| --- | --- | --- |
| Conductividad `k` | 2 | W/(m·K) |
| Densidad `rho` | 2400 | kg/m³ |
| Calor específico `cp` | 900 | J/(kg·K) |

Fuente volumétrica uniforme en el espacio, en W/m³:

```
Q(t) = Qmax/2 · [1 − cos(2πt/τ)]   para 0 ≤ t ≤ τ,      Q(t) = 0 fuera de ese intervalo
```

con **τ = 2 días = 172800 s** y horizonte de evaluación de **48 h** (`TIEMPO_FINAL_S`). El objetivo es el mayor
`Qmax ≥ 0` para que la temperatura del centro, (2, 1) m, no supere **70 °C** en ningún instante del horizonte.

**Método.** Como la temperatura de borde es igual a la inicial, el problema es lineal en `Qmax`:
`T_centro(t) − 20 = Qmax·g(t)`. Se simula una fuente de referencia positiva y se escala:
`Qmax_límite = Qref·(70 − 20) / max_t[T_centro_ref(t) − 20]`. Luego se **verifica con una nueva corrida**.
El máximo se busca en **todos los pasos** del solver (no solo en los fotogramas ni al final).

> **`Qmax` es una aproximación numérica.** Depende de las propiedades del material (`k`, `rho`, `cp`), de la
> malla, del paso de tiempo (Euler explícito) y de la geometría supuesta; no es una garantía exacta del límite
> continuo. El script reporta la sensibilidad a la malla y a `dt` por separado.

## El solver de esta carpeta

`solver_ejercicio.py` es una **copia independiente** de `src/mcplan/thermal/diffusion_nd.py` (commit `8372d4d`)
con ampliaciones solo para este ejercicio (incluida la captura de estados `instantes_captura_s`, que el solver
principal publicado aún no tiene): `fuente_w_m3` puede ser una función del tiempo `f(t [s]) → W/m³`
(evaluada al inicio de cada paso, `t = n·dt`), la energía de la fuente se acumula con lo realmente aplicado,
y `puntos_registro_m` registra la temperatura de puntos físicos en todos los pasos. La copia no importa
`mcplan` y no se mantiene sincronizada con el solver principal; los cambios a uno no afectan al otro.

## Cómo ejecutar (Windows PowerShell o VS Code)

Requisitos: el entorno del proyecto (`.venv` con `numpy` y `matplotlib`; Pillow viene con Matplotlib).

- **Con ▶️ en VS Code:** selecciona el intérprete del proyecto (*Python: Select Interpreter* → `.venv`), abre
  `ejercicio_fuente_variable_2d.py` y pulsa *Ejecutar archivo de Python*. Funciona desde cualquier directorio y no
  requiere `PYTHONPATH`.
- **Desde la terminal** (con el entorno activo, desde la raíz del repositorio):

```powershell
python .\ejercicios_locales\fuente_variable_2d\ejercicio_fuente_variable_2d.py
```

Los parámetros editables están al inicio del script: ventanas y animación (`MOSTRAR_GRAFICOS`,
`MOSTRAR_ANIMACION`, `GUARDAR_ANIMACION`, `NUM_FOTOGRAMAS`, `FPS_ANIMACION`), geometría y malla, propiedades,
temperaturas, `DURACION_FUENTE_S` (= τ, duración de la fuente), `TIEMPO_FINAL_S` (horizonte) y
`QMAX_REFERENCIA_W_M3`. `FPS_ANIMACION` solo controla la velocidad de reproducción. Con
`MOSTRAR_GRAFICOS = False` y `MOSTRAR_ANIMACION = False` se ejecuta sin ventanas (backend Agg) y aun así guarda
las salidas.

## Pruebas locales

```powershell
python -m unittest discover -s .\ejercicios_locales\fuente_variable_2d -v
```

Cubren fuente constante (escalar y función), escalamiento lineal con `Qmax`, balance de energía con fuente
temporal, registro del centro en todos los pasos y comparación del pulso adiabático con
`ΔT(t) = Qmax/(2ρcp)·[t − τ/(2π)·sin(2πt/τ)]` (y `Qmax·τ/(2ρcp)` para `t ≥ τ`), con una tolerancia que incluye el
error de integración temporal. No forman parte de la batería de la entrega (`tests/`).

## Salidas y cómo regenerarlas

Cada ejecución del script escribe en `resultados/` (carpeta ignorada por Git, no versionada):

| Archivo | Contenido |
| --- | --- |
| `resumen.txt` | `Qmax` límite, máximo e instante, centro inicial y final, propiedades, malla, `dt` y límite estable, balance de energía y sensibilidad |
| `centro_cada_paso.csv` | tiempo, `Q` y temperatura del centro en cada paso (`paso, tiempo_s, tiempo_h, Q_inicio_de_paso_W_m3, T_centro_C`) |
| `grafico_fuente_Q_t.png` | `Q(t)` y temperatura del centro con el límite |
| `evolucion_fuente_variable_2d.gif` | mapa térmico y temperatura del centro vs. tiempo, sincronizados (requiere `GUARDAR_ANIMACION = True`) |

Para regenerarlos basta volver a ejecutar el script. `EJECUTAR_SENSIBILIDAD = False` omite las corridas de
sensibilidad (más rápido).
