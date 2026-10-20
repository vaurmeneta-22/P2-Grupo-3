# Reporte semana 1 — Desarrollo del código

**Proyecto:** Proyecto 2 — Planificación de vertidas de hormigón masivo  
**Equipo:** Grupo 3  
**Integrantes:** Monserrat Cubillos Tome, Josefina Ignacia Muro Teuber y Vicente Andrés Urmeneta Campos  
**Repositorio:** `P2-Grupo-3`  
**Fecha de actualización de este registro:** 20 de octubre de 2026

## 1. Propósito y alcance

Durante esta etapa se preparó una primera herramienta numérica para estudiar la difusión del calor en una dimensión. El trabajo permite practicar la ecuación de conducción, comparar condiciones de borde y verificar que el cálculo explícito se mantenga estable.

El programa es una base de aprendizaje y verificación numérica. **No representa todavía la fundación de hormigón masivo del proyecto**, ni debe interpretarse como una predicción de temperatura de una obra real. El modelo actual no incluye hidratación del cemento, fuentes de calor, geometría BIM/IFC, suelo, moldajes, clima del sitio ni incertidumbre.

## 2. Estructura del código relacionado

| Ruta | Responsabilidad |
| --- | --- |
| `src/mcplan/thermal/diffusion_1d.py` | Implementa el solver FTCS, el cálculo de difusividad, las condiciones de borde y las comprobaciones de parámetros y estabilidad. |
| `scripts/demo_difusion_1d.py` | Define un caso editable, genera el perfil inicial, ejecuta simulaciones, grafica perfiles y exporta datos. |
| `scripts/check_ifc.py` | Inspecciona el IFC4 de la fundación, sus vertidas, propiedades, dimensiones y errores de geometría. |
| `models/ifc/Fundacion_1_con_propiedades.ifc` | Modelo BIM de cinco vertidas IFC con propiedades del curso. |
| `models/bim/Fundacion_Grupo3_Presentacion.blend` | Escena Blender de presentación visual vinculada al trabajo BIM. |
| `tests/test_diffusion_1d.py` | Contiene pruebas automáticas de propiedades, condiciones de borde, estabilidad y comportamiento esperado. |
| `requirements.txt` | Fija las dependencias Python de esta etapa: NumPy y Matplotlib. |
| `results/solver_1d/` | Destino que el demo crea para el gráfico y el archivo CSV. Los resultados generados se excluyen de Git por configuración. |

El módulo principal se importa como `mcplan.thermal.diffusion_1d`. El demo añade `src` a la ruta de importación para poder ejecutarse desde la raíz del repositorio.

## 3. Ecuación y propiedades térmicas

El solver implementa difusión térmica sin fuente interna:

\[
\frac{\partial T}{\partial t}
= \alpha\,\frac{\partial^2 T}{\partial x^2},
\qquad
\alpha = \frac{k}{\rho c_p}.
\]

Donde:

- `T` es la temperatura en °C (las diferencias de temperatura son equivalentes en K).
- `t` es el tiempo en segundos.
- `x` es la posición en metros.
- `k` es la conductividad térmica en W/(m K).
- `rho` es la densidad en kg/m³.
- `cp` es el calor específico en J/(kg K).
- `alpha` es la difusividad térmica en m²/s.

`calcular_difusividad_termica(k, rho, cp)` calcula `alpha = k/(rho*cp)` y rechaza propiedades que no sean finitas y positivas. La difusión transporta energía desde zonas de mayor temperatura hacia zonas de menor temperatura; no crea calor por sí sola.

## 4. Discretización espacial y avance temporal

La barra de longitud `L` se representa mediante nodos equiespaciados. El paso espacial es `dx = L/(N - 1)`, donde `N` es el número de nodos. El tiempo se divide en pasos `dt` y se calcula el número de Fourier:

\[
Fo = \frac{\alpha\,dt}{dx^2}.
\]

Para un nodo interior, FTCS actualiza la temperatura usando el estado anterior:

\[
T_i^{n+1}=T_i^n + Fo\left(T_{i+1}^n - 2T_i^n + T_{i-1}^n\right).
\]

El cálculo usa un arreglo NumPy unidimensional. Cada paso se construye desde el arreglo del paso anterior, de modo que no se mezclen valores ya actualizados con valores aún antiguos dentro del mismo paso.

`ResultadoDifusion` devuelve la posición de los nodos, la temperatura al tiempo final, el paso de tiempo realmente utilizado, el número de pasos y el número de Fourier. **El solver devuelve el estado final de una simulación**, no almacena por sí solo la historia completa de todos los pasos.

## 5. Condiciones de borde

La función general `simular_difusion_1d` permite elegir una condición distinta en cada extremo. Los valores asociados a Dirichlet o Neumann deben proporcionarse; en Robin se entregan los parámetros convectivos y la temperatura ambiente.

### 5.1 Dirichlet: temperatura fija

Impone directamente una temperatura en el extremo. Por ejemplo, `T(0,t)=20 °C` mantiene el nodo izquierdo en 20 °C durante la simulación. Se configura con `condicion_*="dirichlet"` y el valor correspondiente en °C.

### 5.2 Neumann: gradiente prescrito

Impone `dT/dx` en K/m. El solver usa una aproximación de punto fantasma para actualizar el nodo de borde. Un gradiente cero equivale a flujo conductivo nulo a través del extremo: el extremo queda aislado. Con ambos extremos aislados, la prueba comprueba que la energía térmica discreta se conserve.

### 5.3 Robin: intercambio convectivo

Representa el intercambio de calor con un ambiente mediante el coeficiente `h` en W/(m² K), la conductividad `k` y una temperatura ambiente. El término convectivo depende de la diferencia entre aire y superficie:

\[
q'' = h\,(T_{aire}-T_{superficie}).
\]

Si el aire está más frío, el borde pierde calor; si está más caliente, el borde puede recibir calor. La actualización de Robin en el código combina el efecto de la conducción desde el nodo vecino con el intercambio convectivo en un volumen de control de borde.

La temperatura ambiente de cada lado puede entregarse como un número constante o como una función del tiempo en segundos. Por ejemplo, una función sinusoidal de media `Tmedia`, amplitud `A` y período `P` tiene la forma:

\[
T_{aire}(t)=T_{media}+A\sin\left(\frac{2\pi t}{P}\right).
\]

Los ambientes izquierdo y derecho son independientes. La amplitud cero deja constante el ambiente de ese lado. Estas funciones solo afectan al solver si la condición del borde respectivo es Robin; con Dirichlet o Neumann no se usan.

## 6. Estabilidad numérica

FTCS explícito requiere un paso temporal limitado. Para difusión 1D con extremos Dirichlet o Neumann, la referencia es `Fo <= 0.5`. Robin añade un término convectivo en el nodo de borde. El código calcula un límite más restrictivo con el número de Biot de malla:

\[
Bi_{malla}=\frac{h\,dx}{k},
\qquad
Fo_{max}=\min\left(0.5,\frac{0.5}{1+Bi_{malla,izq}},
\frac{0.5}{1+Bi_{malla,der}}\right),
\]

considerando únicamente los extremos configurados como Robin. Si el `Fo` resultante supera ese límite, el solver genera un `ValueError` explicando el límite y el valor obtenido. Esta comprobación previene continuar silenciosamente con una solución numérica inestable.

## 7. Qué hace el demo

Al inicio de `scripts/demo_difusion_1d.py` se encuentra el bloque `PARAMETROS EDITABLES`. Entre otros, allí se configuran:

- longitud y número de nodos;
- temperatura inicial del centro y temperatura base;
- tipo de condición de borde y sus valores;
- coeficientes de convección `h`;
- temperatura media, amplitud y período del aire por lado;
- `k`, `rho` y `cp`;
- tiempo final, cantidad e instantes de perfiles y número de Fourier máximo.

`crear_perfil_inicial` crea una curva seno simétrica:

\[
T(x,0)=T_{base}+(T_{centro}-T_{base})\sin\left(\frac{\pi x}{L}\right).
\]

Por tanto, la temperatura base queda en ambos extremos y la temperatura indicada como centro se alcanza a mitad de la barra. Este método **no genera un perfil inicial asimétrico** ni permite definir arbitrariamente la temperatura en cada nodo.

El demo construye `TIEMPOS_GRAFICADOS_S` usando `NUMERO_PERFILES_GRAFICADOS` entre cero y `TIEMPO_FINAL_S`. Para cada instante distinto de cero vuelve a simular desde la condición inicial hasta ese tiempo y toma el estado final; no guarda toda la historia temporal en memoria. El instante cero se grafica como perfil inicial puro. En Dirichlet, esto puede hacer que el perfil de `t=0` muestre nodos extremos distintos de los valores impuestos para `t>0`; el propio script lo indica en un comentario.

Cuando se cumplen las hipótesis del modo seno con temperaturas fijas iguales en ambos extremos, el demo calcula además la solución analítica:

\[
T(x,t)=T_{base}+(T_{centro}-T_{base})\sin\left(\frac{\pi x}{L}\right)
\exp\left[-\alpha\left(\frac{\pi}{L}\right)^2t\right].
\]

En ese caso se puede dibujar la curva analítica final y calcular el error máximo entre la solución numérica y la analítica. En otros casos, esa comparación se omite porque la fórmula no representa las condiciones físicas de la corrida.

## 8. Salidas del demo

Al ejecutarse, el script crea `results/solver_1d/` si no existe y escribe:

1. `perfiles_temperatura_1d.png`: gráfico de temperatura versus posición para los tiempos seleccionados.
2. `temperaturas_1d.csv`: una fila por combinación de tiempo y posición, con las columnas `tiempo_s`, `posicion_m`, `temperatura_numerica_c` y `temperatura_analitica_c`. La columna analítica queda como `NaN` si no aplica la solución de referencia.

El gráfico también se abre en una ventana Matplotlib. Las salidas se pueden regenerar ejecutando nuevamente el demo con la configuración deseada.

## 9. Configuración que está actualmente en el repositorio

La revisión de `scripts/demo_difusion_1d.py` para este reporte encontró:

| Parámetro | Valor actual |
| --- | --- |
| `LONGITUD_M` | 1.0 m |
| `NUMERO_NODOS` | 41 |
| Temperatura inicial centro/base | 0 °C / 0 °C |
| Condición izquierda/derecha | Dirichlet / Dirichlet |
| Temperaturas Dirichlet izquierda/derecha | 20 °C / 20 °C |
| `h` izquierda/derecha | 10 / 10 W/(m² K) (no se usa en esta corrida Dirichlet) |
| Ambiente izquierdo | media 20 °C, amplitud 10 °C, período 12 h (no se usa con Dirichlet) |
| Ambiente derecho | media 20 °C, amplitud 0 °C, período 12 h (no se usa con Dirichlet) |
| `k`, `rho`, `cp` | 2.0 W/(m K), 2000 kg/m³, 10 J/(kg K) |
| `TIEMPO_FINAL_S` | 5000 s |
| Perfiles | 8 |
| `NUMERO_FOURIER_MAXIMO` | 0.45 |

Así, aunque el solver ya admite Robin con ambiente variable, **la configuración por defecto actual no ejecuta esa prueba**. En la corrida actual, los valores de ambiente y los coeficientes `h` no intervienen porque ambos extremos están en Dirichlet. Además, la condición inicial interior es 0 °C y los extremos se imponen a 20 °C al iniciar la integración; hay un salto inicial en los bordes, que debe considerarse al interpretar la curva de `t=0`.

Para realizar una prueba de ambiente variable, revisar y cambiar el bloque de parámetros del demo: seleccionar Robin en el/los extremos que intercambiarán calor, especificar las funciones de aire, elegir `h` y las propiedades térmicas, y comprobar que `Fo` se mantenga bajo el límite. Para que solo varíe el lado izquierdo, usar amplitud cero en el ambiente derecho.

## 10. Pruebas automáticas

Las pruebas en `tests/test_diffusion_1d.py` cubren 14 comportamientos:

1. cálculo de difusividad desde `k`, `rho` y `cp`;
2. comparación con la solución analítica y convergencia al refinar la malla;
3. conservación de un equilibrio uniforme con Dirichlet compatible;
4. imposición de temperaturas fijas en ambos extremos;
5. conservación del promedio trapezoidal con extremos Neumann aislados;
6. permanencia de un estado uniforme con Neumann cero;
7. combinación de Dirichlet y Neumann en una misma barra;
8. rechazo de un paso inestable;
9. rechazo de una temperatura inicial con forma incorrecta;
10. aceptación de una temperatura ambiente que depende del tiempo en Robin;
11. equivalencia entre Robin con `h=0` y Neumann aislado;
12. permanencia del equilibrio cuando la barra y el ambiente están a la misma temperatura;
13. enfriamiento hacia un ambiente más frío;
14. reducción del límite de estabilidad causada por Robin.

Comando usado desde PowerShell en la raíz del repositorio:

```powershell
$env:PYTHONPATH = (Resolve-Path .\src).Path
.\.venv\Scripts\python.exe -m unittest discover -s .\tests -v
```

**Resultado verificado al elaborar/actualizar este reporte:** 14 pruebas ejecutadas, 14 aprobadas (`OK`).

## 11. Cómo preparar el entorno y ejecutar

Desde la raíz del repositorio, en Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
$env:PYTHONPATH = (Resolve-Path .\src).Path
python -m unittest discover -s .\tests -v
python .\scripts\demo_difusion_1d.py
```

Si la política de PowerShell impide activar el entorno, se puede invocar Python directamente:

```powershell
$env:PYTHONPATH = (Resolve-Path .\src).Path
.\.venv\Scripts\python.exe -m unittest discover -s .\tests -v
.\.venv\Scripts\python.exe .\scripts\demo_difusion_1d.py
```

Las dependencias registradas son `numpy==2.5.2` y `matplotlib==3.11.1`.

## 12. Límites conocidos y trabajo pendiente

Esta versión no incorpora todavía:

- fuente volumétrica o calor de hidratación dependiente del tiempo y la temperatura;
- un perfil inicial arbitrario o asimétrico definido por el usuario;
- registro de temperaturas en cada paso temporal del solver (solo devuelve el estado final; el demo vuelve a correr por cada perfil solicitado);
- transferencia de calor en 2D/3D;
- lectura del modelo IFC ni voxelización;
- geometría de la fundación y propiedades por vertida;
- interfaces entre vertidas, suelo, moldajes y aislación;
- datos históricos o pronósticos meteorológicos de Mantos de Hornillos/Ovalle;
- simulación de Monte Carlo, confiabilidad o incertidumbre.

Los valores del demo son ilustrativos y deben identificarse como tales. No deben presentarse como propiedades certificadas de hormigón ni como datos del sitio.

## 13. Relación con el objetivo del Proyecto 2

El solver 1D aporta una base para explicar conducción, balance de energía, discretización, condiciones de borde y estabilidad del método explícito. En las siguientes etapas se debe integrar la geometría y los datos del caso, ampliar la física y verificar cada componente antes de utilizarlo en el análisis de planes de vertidas. El informe final del proyecto deberá dejar trazabilidad entre los parámetros, entradas, código, resultados y SHA de Git evaluado.

## 14. Modelo IFC4 y presentación en Blender/Bonsai

Además del solver didáctico 1D, durante la semana se preparó el primer modelo BIM de la fundación y se organizó una escena visual para presentarlo.

### Modelo IFC

El archivo `models/ifc/Fundacion_1_con_propiedades.ifc` usa esquema IFC4 e incluye cinco elementos `IfcFooting`, uno por vertida. Cada vertida contiene el conjunto `UANDES_HormigonMasivo` con propiedades como `VertidaID`, `Inicio_h`, `TemperaturaColocacion_C`, `TipoHormigon`, `Moldaje`, `Aislacion` y `RetiroMoldaje_h`. La ejecución de `scripts/check_ifc.py` reportó cinco vertidas y cero errores de geometría. El script también imprime las propiedades y las dimensiones calculadas desde la geometría IFC.

Comando desde la raíz del repositorio, con el entorno activo:

```powershell
python .\scripts\check_ifc.py
```

El IFC representa los datos BIM y es la fuente principal para la geometría y propiedades. El solver 1D todavía no lee este modelo ni lo transforma en una malla 3D; esa integración es trabajo pendiente.

### Escena visual guardada en Blender

La presentación está en `models/bim/Fundacion_Grupo3_Presentacion.blend`. La escena conserva cinco vertidas IFC y añade elementos de visualización: textura procedural de hormigón con relieve sutil, tonos ligeramente más oscuros en caras laterales para distinguir los escalones, un plano de arena de 56 × 56 m, iluminación suave y uniforme y vegetación exterior con cactus y arbustos de estepa. La vista queda en perspectiva libre, sin cámara fija, para poder orbitar y navegar por la escena. Estas decisiones visuales no sustituyen ni modifican las propiedades BIM del IFC.

### Qué guardar y cómo mantener ambos archivos

Bonsai y Blender guardan información distinta:

- Para editar las vertidas, sus propiedades o su geometría BIM, se abre el IFC mediante **Bonsai → File → Open IFC Project** y se guarda con **Bonsai → File → Save IFC Project**.
- Para conservar materiales, texturas, luces, arena, vegetación y la organización visual, se guarda la escena con las opciones normales de Blender (**File → Save** o **Save As**) en un `.blend`.
- El `.blend` no se sincroniza automáticamente si el IFC se modifica en otra sesión. Después de actualizar el IFC, hay que cargar la versión vigente mediante Bonsai y actualizar/integrar de nuevo la escena de presentación antes de guardar el `.blend`.
- Conviene conservar ambos: IFC para el modelo BIM y `.blend` para la presentación visual. Las texturas detalladas de render se conservan en Blender, no como textura de imagen dentro del IFC.

Referencias: [guardado y apertura de proyectos IFC en Bonsai](https://docs.bonsaibim.org/studio-unstable/guides/authoring/starting_new_project.html) y [materiales y texturas](https://docs.bonsaibim.org/studio/guides/authoring/advanced_modeling/material_assignment.html).
