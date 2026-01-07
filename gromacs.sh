#!/bin/bash
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

##############################################################################
#                            Global Variables                                #
##############################################################################

CONFIG_FILE="gromacs_config.ini"
LOG_FILE="gromacs_console.log"
BACKUP_DIR="backups"

##############################################################################
#                            Utility Functions                               #
##############################################################################

# Function to log messages with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" &> /dev/null
}

# Function to show error dialog
show_error() {
    dialog --title 'Error' --msgbox "$1" 10 50
    log_message "ERROR: $1"
}

# Function to show info dialog
show_info() {
    dialog --title 'Información' --msgbox "$1" 10 50
    log_message "INFO: $1"
}

# Function to create default configuration file
create_default_config() {
    cat > "$CONFIG_FILE" << 'EOF'
# GROMACS Console Configuration File
# Este archivo contiene todos los parámetros configurables del sistema

[Files]
liganditp = ligand.itp
ligandname = LIG
proteinpdb = protein.pdb
ligandpdb = ligand.pdb

[ForceField]
watermodel = tip3p
forcefield = charmm27

[Box]
distbox = 2
formabox = cubic

[Ions]
pname = NA
nname = CL

[Terminal]
terminal_bg = white
terminal_fg = black
terminal_font = Monospace
terminal_fontsize = 12
terminal_geometry = 92x36

[Paths]
backup_dir = backups
log_file = gromacs_console.log

[Options]
auto_backup = yes
confirm_critical = yes
EOF
    log_message "Created default configuration file: $CONFIG_FILE"
}

# Function to load configuration
load_config() {
    if [ ! -f "$CONFIG_FILE" ]; then
        log_message "Configuration file not found, creating default..."
        create_default_config
        show_info "Se ha creado el archivo de configuración: $CONFIG_FILE\nPuedes editarlo desde el menú principal."
    fi
    
    # Load configuration values
    if [ -f "$CONFIG_FILE" ]; then
        # Parse INI file (simple approach)
        while IFS='=' read -r key value; do
            # Skip comments and empty lines
            [[ $key =~ ^#.*$ ]] && continue
            [[ $key =~ ^\[.*\]$ ]] && continue
            [[ -z $key ]] && continue
            
            # Trim whitespace using bash parameter expansion
            key="${key#"${key%%[![:space:]]*}"}"
            key="${key%"${key##*[![:space:]]}"}"
            value="${value#"${value%%[![:space:]]*}"}"
            value="${value%"${value##*[![:space:]]}"}"
            
            # Export as environment variable
            export "$key=$value"
        done < "$CONFIG_FILE"
        log_message "Configuration loaded successfully"
    fi
}

##############################################################################
#                 Checking si están las cosas instaladas                     #
##############################################################################

# dialog is a utility installed by default on all major Linux distributions.
# But it is good to check availability of dialog utility on your Linux box.

if ! command_exists dialog; then
    echo "Dialog utility is not available, Install it"
    exit 1
fi

if ! command_exists gmx; then
    echo "Gromacs utility is not available, Install it"
    exit 1
fi

if ! command_exists mc; then
    echo "Midnight Commander utility is not available, Install it"
    exit 1
fi

if ! command_exists xmgrace; then
    echo "Grace utility is not available, Install it"
    exit 1
fi

# Load configuration at startup
load_config

# Create MDP template files if they don't exist
create_mdp_files() {
    if [ ! -e "em.mdp" ]; then
        cat > em.mdp << 'EOF'
integrator	= steep
emtol		= 1000.0
emstep      = 0.01
nsteps		= 50000
energygrps	= system

nstlist		    = 20
cutoff-scheme   = Verlet

ns-type		    = grid
rlist		    = 1.0
coulombtype	    = PME
rcoulomb	    = 1.0
rvdw		    = 1.0
pbc             = xyz
EOF
        log_message "Created em.mdp template"
    fi

    if [ ! -e "nvt.mdp" ]; then
        cat > nvt.mdp << 'EOF'
define      = -DPOSRES

integrator  = md
nsteps      = 50000
dt          = 0.002

nstxout     = 5000
nstvout     = 5000
nstenergy   = 5000
nstlog      = 5000
energygrps  = Protein LIG

continuation    = no
constraint-algorithm = lincs
constraints     = all-bonds
lincs-iter      = 1
lincs-order     = 1

cutoff-scheme   = Verlet
ns-type         = grid
nstlist         = 10
rcoulomb        = 1.0
rvdw            = 1.0

coulombtype     = PME
pme-order       = 4
fourierspacing  = 0.135

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1   0.1
ref-t       = 310   310
pcoupl      = no
pbc         = xyz
DispCorr    = EnerPres
gen-vel     = yes
gen-temp    = 310
gen-seed    = -1
EOF
        log_message "Created nvt.mdp template"
    fi

    if [ ! -e "npt.mdp" ]; then
        cat > npt.mdp << 'EOF'
define      = -DPOSRES

integrator  = md
nsteps      = 15000
dt          = 0.002

nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
energygrps  = Protein LIG

continuation    = yes
constraint_algorithm = lincs
constraints     = all-bonds
lincs-iter      = 1
lincs-order     = 2

cutoff-scheme   = Verlet
ns-type         = grid
nstlist         = 20
rcoulomb        = 1.0
rvdw            = 1.0

coulombtype     = PME
pme-order       = 4
fourierspacing  = 0.135

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1   0.1
ref-t       = 310   310

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.58e-5
refcoord-scaling    = com
pbc         = xyz
DispCorr    = EnerPres
gen-vel     = no
EOF
        log_message "Created npt.mdp template"
    fi

    if [ ! -e "md.mdp" ]; then
        cat > md.mdp << 'EOF'
integrator  = md
nsteps      = 20000000
dt          = 0.002

nstxout             = 0
nstvout             = 0
nstenergy           = 10000
nstlog              = 10000
nstxout-compressed  = 50000
compressed-x-grps   = System
energygrps          = Protein LIG

continuation    = yes
constraint_algorithm = lincs
constraints     = all-bonds
lincs-iter      = 1
lincs-order     = 2

cutoff-scheme   = Verlet
ns-type         = grid
nstlist         = 25
rcoulomb        = 1.0
rvdw            = 1.0

coulombtype     = PME
pme-order       = 4
fourierspacing  = 0.12

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1   0.1
ref-t       = 310   310

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.58e-5

pbc         = xyz
DispCorr    = EnerPres
gen-vel     = no
EOF
        log_message "Created md.mdp template"
    fi
}

# Initialize MDP files
create_mdp_files

# Create backup directory if it doesn't exist
[ ! -d "$BACKUP_DIR" ] && mkdir -p "$BACKUP_DIR" && log_message "Created backup directory: $BACKUP_DIR"

##############################################################################
#                      Functions definitions                                 #
##############################################################################

# Helper function to get terminal command with config
get_terminal_cmd() {
    echo "xterm -rightbar -bg ${terminal_bg:-white} -fg ${terminal_fg:-black} -fa '${terminal_font:-Monospace}' -fs ${terminal_fontsize:-12} -geometry ${terminal_geometry:-92x36}"
}

###### edit function: the most basic text editor: NANO #########
function edit() {
    $(get_terminal_cmd) -e "nano -wcSr68; bash"
    log_message "Text editor opened"
}

###### config function: edit configuration file #######
function editconfig() {
    $(get_terminal_cmd) -e "nano -wcSr68 $CONFIG_FILE; bash"
    log_message "Configuration file edited"
    # Reload configuration after editing
    load_config
}

######### topols function: it will create the receptor topology ########
function topols() {
    if [ ! -f "$proteinpdb" ]; then
        show_error "Archivo de proteína '$proteinpdb' no encontrado.\nVerifique la configuración."
        return 1
    fi
    
    $(get_terminal_cmd) -e "gmx pdb2gmx -f $proteinpdb -ff $forcefield -water $watermodel -ignh -o protein-complex.pdb; bash"
    log_message "Topology created for protein: $proteinpdb"
}

########## edtop function #########
function edtop() {
    if [ ! -f "topol.top" ]; then
        show_error "Archivo de topología 'topol.top' no encontrado.\nPrimero debe crear la topología."
        return 1
    fi
    
    $(get_terminal_cmd) -e "nano -wcSr68 topol.top; bash"
    log_message "Topology file edited"
}

########### editp function ############
function editp() {
    if [ ! -f "$liganditp" ]; then
        show_error "Archivo de topología del ligando '$liganditp' no encontrado.\nVerifique la configuración."
        return 1
    fi
    
    $(get_terminal_cmd) -e "nano -wcSr68 $liganditp; bash"
    log_message "Ligand topology file edited"
}

###### neutralise function: crea la caja, solvata y neutraliza ####
function neutral() {
    if [ ! -f "protein-complex.pdb" ]; then
        show_error "Archivo 'protein-complex.pdb' no encontrado.\nPrimero debe crear la topología de la proteína."
        return 1
    fi
    
    $(get_terminal_cmd) -e " \
gmx editconf -f protein-complex.pdb -o protein-complex-box.pdb -c -d $distbox -bt $formabox && \
gmx solvate -cs -cp protein-complex-box.pdb -o protein-complex-solv.pdb -p topol.top && \
gmx grompp -f em.mdp -c protein-complex-solv.pdb -p topol.top -o ions.tpr && \
gmx genion -s ions.tpr -o protein-complex-neutral.pdb -p topol.top -pname $pname -nname $nname -neutral \
; bash"
    log_message "Solvation and neutralization completed"
}

######### function para editar .mdp ##########
function edem() {
    $(get_terminal_cmd) -e "nano -wcSr68 em.mdp; bash"
    log_message "EM parameters edited"
}

function ednvt() {
    $(get_terminal_cmd) -e "nano -wcSr68 nvt.mdp; bash"
    log_message "NVT parameters edited"
}

function ednpt() {
    $(get_terminal_cmd) -e "nano -wcSr68 npt.mdp; bash"
    log_message "NPT parameters edited"
}

function edmd() {
    $(get_terminal_cmd) -e "nano -wcSr68 md.mdp; bash"
    log_message "MD parameters edited"
}

########### corridas function ############

# EM
function emrun() {
    if [ ! -f "protein-complex-neutral.pdb" ]; then
        show_error "Archivo 'protein-complex-neutral.pdb' no encontrado.\nPrimero debe neutralizar el sistema."
        return 1
    fi
    
    $(get_terminal_cmd) -e " \
gmx grompp -f em.mdp -c protein-complex-neutral.pdb -p topol.top -o em.tpr && \
gmx mdrun -v -deffnm em ; bash"
    log_message "Energy minimization completed"
}

# NVT
function nvtrun() {
    if [ ! -f "em.gro" ]; then
        show_error "Archivo 'em.gro' no encontrado.\nPrimero debe ejecutar la minimización de energía."
        return 1
    fi
    
    $(get_terminal_cmd) -e " \
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr && \
gmx mdrun -v -deffnm nvt ; bash"
    log_message "NVT equilibration completed"
}

# NPT
function nptrun() {
    if [ ! -f "nvt.gro" ] || [ ! -f "nvt.cpt" ]; then
        show_error "Archivos 'nvt.gro' o 'nvt.cpt' no encontrados.\nPrimero debe ejecutar el equilibrado NVT."
        return 1
    fi
    
    $(get_terminal_cmd) -e " \
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr && \
gmx mdrun -v -deffnm npt ; bash"
    log_message "NPT equilibration completed"
}

# MD
function mdrun() {
    if [ ! -f "npt.gro" ] || [ ! -f "npt.cpt" ]; then
        show_error "Archivos 'npt.gro' o 'npt.cpt' no encontrados.\nPrimero debe ejecutar el equilibrado NPT."
        return 1
    fi
    
    $(get_terminal_cmd) -e " \
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr && \
gmx mdrun -v -deffnm md ; bash"
    log_message "MD simulation completed"
}

########### New utility functions ############

# Backup function
function backup() {
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    local backup_name="backup_$timestamp"
    local backup_path="$BACKUP_DIR/$backup_name"
    
    mkdir -p "$backup_path"
    
    # Backup important files with error tracking
    local backed_up=0
    local failed=0
    for ext in pdb gro top itp tpr cpt mdp edr log; do
        while IFS= read -r -d '' file; do
            if cp "$file" "$backup_path/"; then
                ((backed_up++))
            else
                ((failed++))
                log_message "WARNING: Failed to backup $file"
            fi
        done < <(find . -maxdepth 1 -name "*.$ext" -print0 2>/dev/null)
    done
    
    if [ -f "$CONFIG_FILE" ]; then
        if cp "$CONFIG_FILE" "$backup_path/"; then
            ((backed_up++))
        else
            ((failed++))
            log_message "WARNING: Failed to backup $CONFIG_FILE"
        fi
    fi
    
    local message="Backup creado en:\n$backup_path\n\nArchivos respaldados: $backed_up"
    if [ $failed -gt 0 ]; then
        message="$message\nFallos: $failed (ver log)"
    fi
    show_info "$message"
    log_message "Backup created: $backup_path ($backed_up files, $failed failed)"
}

# View logs function
function viewlogs() {
    if [ ! -f "$LOG_FILE" ]; then
        show_error "Archivo de log no encontrado."
        return 1
    fi
    
    $(get_terminal_cmd) -e "less $LOG_FILE; bash"
    log_message "Log file viewed"
}

# Cleanup function
function cleanup() {
    dialog --title 'Limpieza de archivos' --yesno 'Esta acción eliminará archivos temporales:\n- Archivos de trayectoria (*.trr, *.xtc grandes)\n- Archivos de energía (*.edr)\n- Archivos de respaldo (#*)\n- mdout.mdp\n\n¿Desea continuar?' 15 60
    
    if [ $? -eq 0 ]; then
        rm -f mdout.mdp 2>/dev/null
        rm -f \#* 2>/dev/null
        
        local count=0
        # Use nullglob to handle case where no files match
        shopt -s nullglob
        for file in *.trr *.edr *.xtc; do
            # Skip compressed trajectory if it's large (probably the production run)
            if [[ "$file" == *.xtc ]] && [[ -f "$file" ]]; then
                # Use portable method to get file size
                local size=$(wc -c < "$file" 2>/dev/null || echo "0")
                # Skip if larger than 100MB (likely production trajectory)
                if [ "$size" -gt 104857600 ]; then
                    continue
                fi
            fi
            rm -f "$file" && ((count++))
        done
        shopt -u nullglob
        
        show_info "Limpieza completada.\n$count archivos eliminados."
        log_message "Cleanup completed: $count files removed"
    fi
}

# Analysis menu function
function analysis() {
    local choice
    dialog --clear --backtitle "Análisis y Visualización" --title "Menú de Análisis" \
        --menu "Seleccione una opción:" 20 70 10 \
        "ENERGY" "Analizar energía (gmx energy)" \
        "RMSD"   "Calcular RMSD" \
        "RMSF"   "Calcular RMSF" \
        "GYRATE" "Radio de giro" \
        "HBOND"  "Puentes de hidrógeno" \
        "BACK"   "Volver al menú principal" 2> menuchoices.$$
    
    choice=$(cat menuchoices.$$)
    
    case $choice in
        ENERGY)
            if [ -f "em.edr" ] || [ -f "nvt.edr" ] || [ -f "npt.edr" ] || [ -f "md.edr" ]; then
                $(get_terminal_cmd) -e "echo 'Archivos .edr disponibles:' && ls -1 *.edr 2>/dev/null && echo -e '\nEjemplo: gmx energy -f em.edr -o energy.xvg' && bash"
            else
                show_error "No se encontraron archivos .edr"
            fi
            ;;
        RMSD)
            if [ -f "md.tpr" ]; then
                $(get_terminal_cmd) -e "gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns; xmgrace rmsd.xvg & bash"
                log_message "RMSD analysis performed"
            else
                show_error "Archivo md.tpr no encontrado"
            fi
            ;;
        RMSF)
            if [ -f "md.tpr" ]; then
                $(get_terminal_cmd) -e "gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res; xmgrace rmsf.xvg & bash"
                log_message "RMSF analysis performed"
            else
                show_error "Archivo md.tpr no encontrado"
            fi
            ;;
        GYRATE)
            if [ -f "md.tpr" ]; then
                $(get_terminal_cmd) -e "gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg; xmgrace gyrate.xvg & bash"
                log_message "Gyration analysis performed"
            else
                show_error "Archivo md.tpr no encontrado"
            fi
            ;;
        HBOND)
            if [ -f "md.tpr" ]; then
                $(get_terminal_cmd) -e "gmx hbond -s md.tpr -f md.xtc -o hbond.xvg -num; xmgrace hbond.xvg & bash"
                log_message "H-bond analysis performed"
            else
                show_error "Archivo md.tpr no encontrado"
            fi
            ;;
        BACK)
            # Return to main menu
            return 0
            ;;
    esac
}

# File manager function
function filemanager() {
    if command_exists mc; then
        mc
        log_message "File manager opened"
    else
        show_error "Midnight Commander no está instalado"
    fi
}

# Help function
function showhelp() {
    dialog --title 'Ayuda - GROMACS Console' --msgbox '\
GROMACS Console - Sistema de Interfaz para Simulaciones MD

FLUJO DE TRABAJO TÍPICO:
1. CONFIG: Editar configuración (archivos, parámetros)
2. TOPOLS: Crear topología de la proteína
3. EDTOP/EDITP: Editar topologías si es necesario
4. NEUTRAL: Solvatar y neutralizar el sistema
5. EDEM -> EM: Configurar y ejecutar minimización
6. EDNVT -> NVT: Configurar y ejecutar equilibrado NVT
7. EDNPT -> NPT: Configurar y ejecutar equilibrado NPT
8. EDMD -> MD: Configurar y ejecutar dinámica molecular
9. ANALYSIS: Analizar resultados

FUNCIONES ADICIONALES:
- BACKUP: Crear copia de seguridad de archivos
- CLEANUP: Limpiar archivos temporales
- LOGS: Ver registro de operaciones
- FILES: Explorador de archivos (Midnight Commander)

CONFIGURACIÓN:
El archivo gromacs_config.ini contiene todos los parámetros
configurables del sistema. Se crea automáticamente si no existe.

Para más información: README.md' 30 75
    log_message "Help displayed"
}


##############################################################################
#                          MENU DE LA APLICACION                             #
##############################################################################

trap 'rm -f menuchoices.*' EXIT     # borrar el archivo temporal

# Display welcome message on first run
if [ ! -f ".gromacs_console_initialized" ]; then
    dialog --title 'Bienvenido a GROMACS Console' --msgbox '\
¡Bienvenido a GROMACS Console!\n\n\
Este sistema proporciona una interfaz simplificada\n\
para realizar simulaciones de dinámica molecular.\n\n\
Se ha creado un archivo de configuración inicial.\n\
Puede editarlo desde CONFIG en el menú principal.\n\n\
Presione OK para continuar...' 15 60
    touch .gromacs_console_initialized
    log_message "GROMACS Console initialized"
fi

while :
do
    # Dialog utility to display options list
    dialog --clear --backtitle "GROMACS Console v2.0 - Sistema Optimizado de Simulación MD" \
    --title "Menú Principal" \
    --menu "Use las flechas para navegar, Enter para seleccionar:" 25 80 18 \
    "CONFIG"   "Editar archivo de configuración (gromacs_config.ini)" \
    "TOPOLS"   "Crear topología de la proteína (topol.top)" \
    "EDTOP"    "Editar topología de la proteína" \
    "EDITP"    "Editar topología del ligando (.itp)" \
    "NEUTRAL"  "Solvatación y neutralización del sistema" \
    "EDEM"     "Editar parámetros de minimización (em.mdp)" \
    "EM"       "Ejecutar minimización de energía" \
    "EDNVT"    "Editar parámetros NVT (nvt.mdp)" \
    "NVT"      "Ejecutar equilibrado NVT (T constante)" \
    "EDNPT"    "Editar parámetros NPT (npt.mdp)" \
    "NPT"      "Ejecutar equilibrado NPT (P constante)" \
    "EDMD"     "Editar parámetros MD (md.mdp)" \
    "MD"       "Ejecutar dinámica molecular" \
    "ANALYSIS" "Análisis y visualización de resultados" \
    "BACKUP"   "Crear copia de seguridad" \
    "CLEANUP"  "Limpiar archivos temporales" \
    "FILES"    "Explorador de archivos (Midnight Commander)" \
    "LOGS"     "Ver registro de operaciones" \
    "EDIT"     "Editor de texto básico (nano)" \
    "HELP"     "Ayuda y documentación" \
    "EXIT"     "Salir del programa" 2> menuchoices.$$

    retopt=$?
    choice=$(cat menuchoices.$$ 2>/dev/null)

    case $retopt in
        0)
            case $choice in
                CONFIG)   editconfig ;;
                EDIT)     edit ;;
                TOPOLS)   topols ;;
                EDTOP)    edtop ;;
                EDITP)    editp ;;
                NEUTRAL)  neutral ;;
                EDEM)     edem ;;
                EDNVT)    ednvt ;;
                EDNPT)    ednpt ;;
                EDMD)     edmd ;;
                EM)       emrun ;;
                NVT)      nvtrun ;;
                NPT)      nptrun ;;
                MD)       mdrun ;;
                ANALYSIS) analysis ;;
                BACKUP)   backup ;;
                CLEANUP)  cleanup ;;
                FILES)    filemanager ;;
                LOGS)     viewlogs ;;
                HELP)     showhelp ;;
                EXIT)
                    log_message "GROMACS Console session ended"
                    clear
                    echo "¡Gracias por usar GROMACS Console!"
                    exit 0
                    ;;
            esac
            ;;
        *)
            log_message "GROMACS Console session cancelled"
            clear
            exit
            ;;
    esac
done
