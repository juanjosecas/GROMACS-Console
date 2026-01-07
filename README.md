# GROMACS Console

## Antecedentes y objetivo: 

Las simulaciones de dinámica molecular (MDS) utilizando GROMACS se encuentran entre los experimentos computacionales comúnmente utilizados en el área de la biología molecular y el descubrimiento de fármacos.

## Materiales y métodos:

La herramienta se construyó utilizando la programación bash shell, y la interfaz de texto de usuario se construyó utilizando el software 'dialog'.

## Resultados y conclusiones:

Esta herramienta ofrece un marco simple, semiautomatizado y relativamente rápido para lo que antes era una tarea compleja, manual, lenta y propensa a errores, presentando un método útil para bioquímicos y químicos sintéticos sin experiencia previa en la interfaz de línea de comandos.

## Características principales (v2.0)

### Sistema de configuración externo
- Archivo de configuración `gromacs_config.ini` con todos los parámetros configurables
- Creación automática del archivo de configuración si no existe
- Configuración de archivos de entrada, campos de fuerza, parámetros de caja, iones, terminal y más

### Funciones optimizadas
- **Validación de archivos**: Verificación de existencia de archivos antes de ejecutar comandos
- **Mensajes de error mejorados**: Indicaciones claras sobre qué archivos faltan y qué hacer
- **Registro de operaciones**: Todas las operaciones se registran en `gromacs_console.log`
- **Funciones auxiliares**: Utilidades para logging, validación y manejo de errores

### Nuevas funciones añadidas
1. **BACKUP**: Crea copias de seguridad automáticas de archivos importantes (.pdb, .gro, .top, .itp, .mdp, etc.)
2. **CLEANUP**: Limpieza inteligente de archivos temporales (con confirmación)
3. **ANALYSIS**: Menú de análisis con opciones para:
   - Análisis de energía
   - Cálculo de RMSD
   - Cálculo de RMSF
   - Radio de giro
   - Puentes de hidrógeno
4. **FILES**: Acceso rápido al explorador de archivos (Midnight Commander)
5. **LOGS**: Visualización del registro de operaciones
6. **HELP**: Sistema de ayuda integrado con flujo de trabajo
7. **CONFIG**: Editor del archivo de configuración externo

### Interfaz mejorada
- Título descriptivo con número de versión
- Descripciones más claras de cada opción
- Mensaje de bienvenida en el primer uso
- Mejor organización del menú
- Mensajes de confirmación para operaciones críticas

## Uso

### Requisitos previos
- GROMACS (gmx)
- Dialog
- Midnight Commander (mc)
- XMGrace (xmgrace)
- XTerm

### Instalación y ejecución
```bash
# Hacer el script ejecutable
chmod +x gromacs.sh

# Ejecutar la consola
./gromacs.sh
```

### Primera ejecución
En la primera ejecución, el script:
1. Verificará las dependencias necesarias
2. Creará el archivo de configuración `gromacs_config.ini`
3. Creará los archivos de plantilla MDP si no existen
4. Mostrará un mensaje de bienvenida
5. Creará el directorio de backups

### Flujo de trabajo típico
1. **CONFIG**: Editar configuración (nombres de archivos, parámetros)
2. **TOPOLS**: Crear topología de la proteína
3. **EDTOP/EDITP**: Editar topologías si es necesario
4. **NEUTRAL**: Solvatar y neutralizar el sistema
5. **EDEM** → **EM**: Configurar y ejecutar minimización de energía
6. **EDNVT** → **NVT**: Configurar y ejecutar equilibrado NVT
7. **EDNPT** → **NPT**: Configurar y ejecutar equilibrado NPT
8. **EDMD** → **MD**: Configurar y ejecutar dinámica molecular
9. **ANALYSIS**: Analizar resultados

### Funciones de utilidad
- **BACKUP**: Recomendado antes de cambios importantes
- **CLEANUP**: Limpiar archivos grandes temporales
- **FILES**: Navegar y gestionar archivos fácilmente
- **LOGS**: Revisar el historial de operaciones

## Estructura de archivos

```
.
├── gromacs.sh                  # Script principal
├── gromacs_config.ini          # Configuración (creado automáticamente)
├── *.mdp                       # Archivos de parámetros MD (creados automáticamente)
├── gromacs_console.log         # Registro de operaciones
└── backups/                    # Directorio de copias de seguridad
    └── backup_YYYYMMDD_HHMMSS/ # Backups con timestamp
```

## Configuración

El archivo `gromacs_config.ini` contiene las siguientes secciones:

- **[Files]**: Nombres de archivos de entrada (proteína, ligando, etc.)
- **[ForceField]**: Campo de fuerza y modelo de agua
- **[Box]**: Parámetros de la caja de simulación
- **[Ions]**: Nombres de iones para neutralización
- **[Terminal]**: Configuración de la terminal (colores, fuente, tamaño)
- **[Paths]**: Rutas para backups y logs
- **[Options]**: Opciones adicionales

## Mejoras técnicas

### Optimizaciones del código
- Uso de funciones auxiliares para reducir duplicación de código
- Funciones de validación centralizadas
- Sistema de logging consistente
- Manejo robusto de errores
- Mejor organización del código con comentarios

### Mantenibilidad
- Configuración separada del código
- Plantillas de archivos MDP más limpias
- Código modular y reutilizable
- Nombres de funciones descriptivos

## Licencia

GNU General Public License v2.0 or later

## Contribuciones

Este proyecto está abierto a contribuciones. Por favor, siéntase libre de reportar problemas o sugerir mejoras.
