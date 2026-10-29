# Ejercicio 3D: fuente variable en el tiempo y cuatro tipos de borde

Ejercicio de clases **independiente** de E04/E05 y de los ejercicios 2D: no modifica ni depende de `src/mcplan/`,
`tests/`, los demos, `scripts/verificar_E05.py` ni `ejercicios_locales/fuente_variable_2d/`. Usa su propia copia del
solver (`solver_ejercicio.py`). Es un ejercicio de clases, no una entrega.

## Problema físico

Caja de **3 m (x) × 3 m (y) × 0,60 m (z)**, temperatura inicial uniforme de **20 °C** en todos los escenarios,
material homogéneo de propiedades constantes (los del ejercicio 2D):

| Propiedad | Valor | Unidad |
| --- | --- | --- |
| Conductividad `k` | 2 | W/(m·K) |
| Densidad `rho` | 2400 | kg/m³ |
| Calor específico `cp` | 900 | J/(kg·K) |

Fuente uniforme en el espacio: `Q(t) = Qmax/2·[1 − cos(2πt/τ)]` para `0 ≤ t ≤ τ`, y `0` después, con
**τ = 172800 s** y horizonte de **48 h** (ambos editables e independientes).

**Qmax = 1282,5659131507364 W/m³**, el valor del ejercicio Dirichlet 2D. Procedencia: el resumen 2D publicaba
1282.565913; el valor completo se recuperó de `fuente_variable_2d/resultados/centro_cada_paso.csv`, despejando
`Q(t_n)` en 312 pasos (dispersión 5e-14). Ese CSV no se versiona (se regenera ejecutando el ejercicio 2D); el valor
queda fijado en `QMAX_W_M3` y una prueba comprueba la coincidencia con el CSV cuando existe. **No se recalcula ni se
optimiza** para esta geometría; superar 70 °C es un resultado a informar.

## Escenarios (cada uno es una simulación aparte)

1. Dirichlet 20 °C en las seis caras.
2. Dirichlet 30 °C en las seis caras (el interior parte a 20 °C; el borde vale 30 °C desde `t = 0`).
3. Neumann con flujo cero (aislamiento).
4. Robin en las seis caras, ambiente a 20 °C, para una lista editable de `beta` (por defecto 0,1; 1; 10 1/m).

**Convención Robin.** `beta = h/k` [1/m] y `n_ext` es la normal **exterior**:
`−k·∂T/∂n_ext = h·(T_s − T_amb)`, es decir `∂T/∂n_ext = −beta·(T_s − T_amb)`. Si el cuerpo está más caliente que el
ambiente, el calor sale (se enfría); si está más frío, entra. El pizarrón muestra el signo positivo: corresponde a
derivar respecto de la normal **interior** (`∂T/∂n_int = +beta·(T_s − T_amb)`), que es lo mismo. El solver usa
`robin(h, T_ambiente)`, por eso `h = k·beta` (con `k = 2`: `h = 0,2; 2; 20 W/(m²·K)`). **Los `beta` son valores
exploratorios, no datos entregados por el profesor.**

## Solver y cálculo

`solver_ejercicio.py` es una copia independiente del solver del ejercicio 2D (que a su vez copia el del repositorio,
commit `8372d4d`), con ampliaciones para este ejercicio: fuente como función del tiempo evaluada al inicio de cada
paso, energía de la fuente realmente aplicada, temperatura interpolada en puntos físicos en todos los pasos y
(`registrar_metricas`) máximo, mínimo y media del dominio en todos los pasos. Opera en 3D (celdas de volumen finito;
caras, aristas y esquinas suman las contribuciones de sus caras expuestas).

- Malla inicial 60×60×12 (dx = dy = dz = 0,05 m); centro físico (1,5; 1,5; 0,30) m por interpolación.
- **dt común** para comparar escenarios: `FRACCION_PASO_ESTABLE` (0,8) del límite del escenario más exigente (el
  Dirichlet, cuya esquina con tres caras de temperatura fija manda: `rho·cp·dx²/(9k) = 300 s`) y con un número entero
  de pasos por fotograma (por defecto 225,29 s).
- Se guardan campos completos solo en los fotogramas y de ellos solo los tres cortes; las series (centro, máximo,
  mínimo y media del dominio) se registran en todos los pasos.
- Los cortes pasan por el centro: `xy` en `z = 0,30 m`, `xz` en `y = 1,5 m` y `yz` en `x = 1,5 m`. Cuando el plano cae
  entre dos capas de celdas (p. ej. `z = 0,30 m` entre capas a 0,275 y 0,325 m) se interpola linealmente; la
  animación y el resumen indican las coordenadas reales. La escala de color es fija y común a todos los escenarios.

## Cómo ejecutar (Windows PowerShell o VS Code)

Requisitos: el entorno del proyecto (`.venv` con `numpy` y `matplotlib`; Pillow viene con Matplotlib). No hace falta
instalar nada nuevo.

- **Con ▶️ en VS Code:** selecciona el intérprete `.venv`, abre `ejercicio_fuente_variable_3d.py` y pulsa
  *Ejecutar archivo de Python*. Funciona desde cualquier directorio y no requiere `PYTHONPATH`.
- **Desde la terminal** (con el entorno activo, desde la raíz del repositorio):

```powershell
python .\ejercicios_locales\fuente_variable_3d\ejercicio_fuente_variable_3d.py
python -m unittest discover -s .\ejercicios_locales\fuente_variable_3d -v
```

Opciones al inicio del script: `MOSTRAR_GRAFICOS`, `MOSTRAR_ANIMACION`, `GUARDAR_ANIMACION`,
`ESCENARIOS_ANIMADOS`, `NUM_FOTOGRAMAS`, `FPS_ANIMACION` (solo reproducción), `EJECUTAR_SENSIBILIDAD`, geometría y
malla, propiedades, `QMAX_W_M3`, `DURACION_FUENTE_S` (τ), `TIEMPO_FINAL_S`, `BETAS_ROBIN_1_M`. Con las dos opciones de
mostrar en `False` se ejecuta sin ventanas (Agg) y aun así guarda las salidas. Hay una sola llamada final a
`plt.show()`. La ejecución completa (con sensibilidad y los GIF) tarda unos minutos.

## Salidas (`resultados/`, no se versionan)

- `comparacion/`: `resumen.txt` (parámetros, tabla comparativa, balance de energía, verificaciones y sensibilidad),
  `tabla_comparativa.csv` y `grafico_comparativo.png` (Q(t) común y temperatura central de todos los escenarios con la
  línea de 70 °C).
- `<escenario>/` (`dirichlet_20C`, `dirichlet_30C`, `neumann_cero`, `robin_beta_0p1`, …): `resumen.txt`,
  `centro_cada_paso.csv` (centro y máximo, mínimo y media del dominio en cada paso) y `evolucion_<escenario>.gif`
  (tres cortes sincronizados con la curva central).

## Verificaciones y pruebas

El script ejecuta, en malla gruesa, Robin con `h = 0` frente al aislamiento, el enfriamiento de un cuerpo caliente
con Robin sin fuente, el equilibrio con Dirichlet igual a la temperatura inicial, el caso adiabático frente a
`T(t) = T0 + Qmax/(2ρcp)·[t − τ/(2π)·sin(2πt/τ)]` (≈ 71,30 °C a las 48 h) y el balance de energía; además informa
simetría, uniformidad en Neumann cero y sensibilidad de malla y de dt (con su costo). Las mismas comprobaciones
están en `test_fuente_variable_3d.py`.

## Limitaciones

- `Qmax` es el del caso 2D Dirichlet; en 3D no está calibrado para esta geometría.
- Resultados numéricos (volúmenes finitos, Euler explícito): dependen de la malla y de `dt`; las diferencias de
  sensibilidad son indicadores, no cotas de error.
- Los `beta` Robin son exploratorios; las propiedades son constantes y homogéneas; no hay hidratación cinética,
  clima ni geometría IFC.
