# Proyecto 2 — Planificación de vertidas de hormigón masivo

Proyecto de Métodos Computacionales · Grupo 3

## Equipo

- Monserrat Cubillos Tome
- Josefina Ignacia Muro Teuber
- Vicente Andres Urmeneta Campos

## Objetivo

Diseñar, justificar y defender un plan de construcción por vertidas para una fundación de hormigón masivo, considerando el comportamiento térmico y la incertidumbre de los materiales, la ejecución y el clima.

El producto final será una decisión constructiva respaldada por simulaciones y un Plan de Control de Temperatura. El horizonte de análisis es de **14 días desde el inicio de la primera vertida**. La fecha hipotética de construcción corresponde a la entrega/defensa final, todavía por definir.

## Sitio asignado

| Campo | Valor |
| --- | --- |
| Grupo | 3 |
| Sitio | Mantos de Hornillos / Ovalle |
| Región | Coquimbo |
| Parque | Camarico – Ecos del Mar |
| Estado indicado en el curso | Ampliación/proyecto |
| Latitud | -31.111500 |
| Longitud | -71.614700 |

Todos los scripts meteorológicos deben usar exactamente estas coordenadas, salvo autorización explícita del profesor. Son el punto meteorológico oficial del caso docente; no representan necesariamente el centro de una fundación real.

Las características climáticas esperadas del sector son hipótesis que deben contrastarse con datos reales.

## Caso base y criterios térmicos

La geometría se basa en la fundación de turbina eólica estudiada por Ju y Lei (2019): dimensiones aproximadas de 18 m × 18 m, espesor de 3,30 m y una parte superior cilíndrica de aproximadamente 7 m de diámetro. Se utilizará la geometría oficial entregada o validada por el profesor. Las dimensiones faltantes y simplificaciones deberán documentarse.

Criterios del proyecto docente:

- Temperatura máxima: límite de **70 °C**.
- Diferencial núcleo–superficie: límite fijo de **20 °C** durante todo el análisis.
- Temperatura superficial evaluada en puntos/celdas representativos aproximadamente **25 mm bajo las superficies exteriores**.

Cumplir estos criterios se interpreta como control del riesgo de agrietamiento térmico, no como una demostración mecánica de ausencia de grietas. Los límites se adoptan para este curso y no deben presentarse como exigencias universales de la especificación ICH.

Las ecuaciones exactas de evaluación y confiabilidad se transcribirán desde el enunciado original: no quedaron visibles en la copia de texto recibida.

## Flujo de trabajo previsto

```text
Bonsai / Blender → IFC → IfcOpenShell → geometría y propiedades
                                            ↓
                                       voxelización
                                            ↓
                                  solver térmico 3D propio
                                            ↓
                          métricas térmicas y confiabilidad
                                            ↓
                       comparación de planes y recomendación
```

El solver será independiente del estado interno de Blender y utilizará una grilla cartesiana regular con balance térmico por celda. Se desarrollará por el equipo; un solver externo no sustituirá el núcleo numérico.

La física mínima incluirá conducción, propiedades térmicas, calor de hidratación, temperatura de colocación, activación temporal e interfaces de vertidas, convección, interacción con suelo, moldajes, aislación, retiro de protecciones y temperatura ambiente variable.

Cada vertida deberá ser identificable en el IFC, con información equivalente a `VertidaID`, `TipoHormigon`, `InicioPlanificado`, `TemperaturaColocacion_C`, `Moldaje`, `Aislacion` y `RetiroMoldaje_h`.

## Estructura del repositorio

```text
P2-Grupo-3/
├── README.md
├── .gitignore
├── configs/                  # Casos, parámetros y planes de vertidas
├── data/
│   ├── raw/                  # Datos originales y sus metadatos
│   ├── weather/
│   │   ├── historical/       # Series meteorológicas históricas
│   │   └── forecasts/        # Snapshots de pronóstico, sin sobrescritura
│   └── processed/            # Datos preparados para simulación
├── models/
│   ├── bim/                  # Modelos fuente CAD/BIM
│   └── ifc/                  # Exportaciones IFC para el solver
├── src/
│   └── mcplan/
│       ├── geometry/         # Lectura IFC y voxelización
│       ├── thermal/          # Solver, hidratación y condiciones de borde
│       ├── weather/          # Adquisición y preparación del clima
│       ├── reliability/      # Muestreo, fallas y agregación estadística
│       └── visualization/    # Campos 3D, curvas y gráficos
├── tests/                    # Verificación numérica y pruebas del software
├── scripts/                  # Utilidades y ejecución reproducible
├── hpc/                      # Jobs Slurm y ejecución por lotes
├── notebooks/                # Exploración y análisis documentado
├── results/                  # Salidas generadas; no se versionan por defecto
└── docs/
    ├── references/           # Referencias y antecedentes permitidos
    ├── decisions/            # Supuestos y decisiones de modelación
    └── deliverables/         # Evidencias y entregas del curso
```

Las carpetas vacías se conservan en Git mediante archivos `.gitkeep`. La estructura es inicial y puede ajustarse a medida que se implementen los módulos.

## Estado actual y entorno

**Avance E01.** El solver 1D FTCS permite seleccionar Dirichlet, Neumann o Robin independientemente en cada extremo y cuenta con pruebas automatizadas; NumPy está fijado en `requirements.txt`. Neumann se expresa como gradiente `dT/dx` en K/m; un gradiente cero representa un extremo aislado. Robin modela convección con temperatura ambiente constante o variable en el tiempo y coeficiente `h`.

También se preparó un modelo IFC4 de la fundación con cinco vertidas `IfcFooting` y propiedades del curso. `scripts/check_ifc.py` inspecciona sus propiedades y dimensiones; la última ejecución registrada encontró cinco vertidas y cero errores de geometría. La escena visual se conserva en `models/bim/Fundacion_Grupo3_Presentacion.blend`. El solver todavía no importa ni voxeliza el IFC; faltan calor de hidratación, geometría 3D y el flujo completo de simulación.

Herramientas previstas: Python, Git/GitHub, Bonsai/Blender, IfcOpenShell, PyVista/VTK o ParaView, Matplotlib y acceso al cluster del curso. Las versiones y los pasos de instalación se documentarán al configurar el entorno.

La interfaz final deberá permitir operaciones equivalentes a las siguientes; estos comandos son objetivos y **aún no están disponibles**:

```text
mcplan simulate plan.yaml
mcplan reliability plan.yaml
mcplan compare plan_A.yaml plan_B.yaml
```

## Datos, incertidumbre y reproducibilidad

- Conservar los datos originales con su procedencia y unidades.
- Guardar cada pronóstico como un snapshot independiente con fuente, modelo meteorológico, fecha/hora de descarga, coordenadas, unidades, horizonte y variables. No reemplazar silenciosamente versiones anteriores.
- Explicitar cómo se completa el horizonte de 14 días cuando el pronóstico no lo cubra y actualizar la información meteorológica para la entrega final.
- Distinguir incertidumbre de materiales, ejecución, clima y, cuando corresponda, del modelo. Justificar las distribuciones utilizadas.
- Reportar la incertidumbre estadística de las estimaciones de confiabilidad.
- Registrar configuración, semillas aleatorias cuando correspondan, datos de entrada y commit de cada corrida.
- Ejecutar las realizaciones independientes de Monte Carlo mediante arrays o lotes en HPC.
- Mantener los resultados voluminosos fuera de Git; documentar su ubicación y los pasos para regenerarlos.

## Entregas y evaluación

Toda evidencia de código se evalúa sobre un SHA explícito. Una vez iniciada la revisión, ese commit queda congelado para la actividad.

- Repositorio: https://github.com/vaurmeneta-22/P2-Grupo-3
- Commit inicial de creación: `b79ed3a31d9e3689b4a6c8e41ee831c102523bb5`
- E00: jueves 15 de octubre de 2026, 10:30; 0 puntos, requisito para actividades posteriores.
- Fecha de entrega/defensa final: pendiente de confirmación.

El commit inicial antecede a esta estructura. Si se entrega una versión posterior, se deberá declarar el SHA que efectivamente contenga esa evidencia.

El Plan de Control de Temperatura final incluirá materiales, hidratación, temperaturas previstas, procedimiento de vertidas, protecciones, monitoreo y sensores, retiro de protecciones, contingencias y confiabilidad. Cada integrante deberá poder explicar el código, las ecuaciones, los datos y las decisiones del equipo.

## Antecedentes pendientes

- Geometría oficial del curso.
- Texto completo de las ecuaciones omitidas en la copia del enunciado.
- Paper de Ju y Lei (2019) y especificación ICH del curso para su revisión detallada.
- Fecha final y requisitos de las próximas entregas.

## Solver térmico 1D — etapa inicial E01

Esta etapa implementa difusión térmica 1D con el esquema explícito FTCS y permite elegir Dirichlet, Neumann o Robin en cada extremo. Para Robin, la temperatura ambiente puede ser constante o una función del tiempo. El demo permite definir promedio, amplitud y período independientemente para cada lado: el izquierdo tiene 20 ± 10 °C con período de 12 h, y el derecho queda constante a 20 °C porque su amplitud está en cero. Con Dirichlet, `VALOR_IZQUIERDO` / `VALOR_DERECHO` son temperaturas fijas en °C; con Neumann son gradientes `dT/dx` en K/m, donde `0.0` modela aislamiento; con Robin se usan las temperaturas ambiente y `COEF_CONVECCION_IZQUIERDO_W_M2_K` / `COEF_CONVECCION_DERECHO_W_M2_K` en W/(m² K). El solver calcula la difusividad desde `k`, `rho` y `cp`, informa el número de Fourier y comprueba la estabilidad considerando también Robin. El caso seno de referencia se compara con su solución analítica cuando corresponden esas hipótesis.

Todavía no incluye calor de hidratación, interfaces, geometría IFC ni clima. Es un caso de verificación numérica inicial, no una simulación de la fundación ni un plan constructivo.

### Preparar entorno y ejecutar pruebas (Windows PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
$env:PYTHONPATH = "src"
python -m unittest discover -s tests -v
```

El archivo `requirements.txt` fija la versión de NumPy usada por esta etapa. Se requiere Python 3.12 o compatible con esa versión.

### Graficar y configurar el ejemplo del solver 1D

Con el entorno virtual activo, ejecuta:

```powershell
python .\scripts\demo_difusion_1d.py
```

Antes de ejecutarlo, abre `scripts/demo_difusion_1d.py` y modifica la sección `PARAMETROS EDITABLES` que aparece al principio. Allí puedes cambiar las temperaturas iniciales del centro y los extremos, las temperaturas fijas de ambos bordes, las propiedades `k`, `rho` y `cp`, la longitud, el número de nodos, el tiempo final y los instantes que aparecerán en el gráfico. El programa calcula `alfa = k/(rho*cp)` automáticamente.

Al ejecutarlo, el script muestra los perfiles de temperatura y guarda `perfiles_temperatura_1d.png` y `temperaturas_1d.csv` en `results/solver_1d/`. La comparación con la solución analítica se muestra cuando se mantienen las condiciones del caso seno de referencia. Este ejemplo es de difusión pura; no representa todavía la fundación ni una mezcla real.

## Modelo BIM/IFC y escena de presentación

El modelo BIM y su presentación visual son archivos relacionados, pero no intercambiables:

| Archivo | Uso |
| --- | --- |
| `models/ifc/Fundacion_1_con_propiedades.ifc` | Modelo IFC4 para el flujo BIM, con cinco vertidas `IfcFooting` y el conjunto de propiedades `UANDES_HormigonMasivo`. Incluye, entre otros, identificador de vertida, hora de inicio, temperatura de colocación, tipo de hormigón, moldaje, aislación y retiro de moldaje. |
| `models/bim/Fundacion_Grupo3_Presentacion.blend` | Escena de Blender para presentación visual: materiales y textura de hormigón, terreno de arena, iluminación, vegetación y vista interactiva. |
| `scripts/check_ifc.py` | Inspección del IFC: esquema, propiedades, dimensiones y errores de geometría. |

La presentación Blender contiene una textura de hormigón de bajo relieve para dejar visibles los niveles, caras laterales ligeramente más oscuras, un plano de arena de 56 × 56 m, iluminación suave y uniforme, y cactus y arbustos en el exterior. No se añadió una cámara fija: la vista permanece en perspectiva libre para orbitar y navegar. Estas adiciones son de visualización y no cambian la geometría IFC de las vertidas.

### Guardado y actualización

- Para cambios BIM, abrir `Fundacion_1_con_propiedades.ifc` con **Bonsai → File → Open IFC Project** y guardar con **Bonsai → File → Save IFC Project**.
- Para conservar la escena de presentación, abrir `Fundacion_Grupo3_Presentacion.blend` con las opciones normales de Blender y guardar con **File → Save** o **Save As**.
- El `.blend` no se actualiza automáticamente cuando cambia el IFC. Si el IFC se modifica en otra sesión, se debe cargar la versión actualizada en Bonsai y volver a integrar la presentación visual antes de guardar el `.blend`.
- Las texturas detalladas de Blender se guardan en el `.blend`; el IFC sigue siendo la fuente principal de la geometría y los datos BIM.

Bonsai separa el guardado IFC del guardado Blender. Ver la [guía de guardado de proyectos IFC](https://docs.bonsaibim.org/studio-unstable/guides/authoring/starting_new_project.html) y la [guía de materiales y texturas](https://docs.bonsaibim.org/studio/guides/authoring/advanced_modeling/material_assignment.html).
