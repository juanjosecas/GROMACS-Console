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
#                 Checking si están las cosas instaladas                     #
##############################################################################

# dialog is a utility installed by default on all major Linux distributions.
# But it is good to check availability of dialog utility on your Linux box.

which dialog &> /dev/null

[ $? -ne 0 ]  && echo "Dialog utility is not available, Install it" && exit 1

which gmx &> /dev/null

[ $? -ne 0 ]  && echo "Gromacs utility is not available, Install it" && exit 1

which mc &> /dev/null

[ $? -ne 0 ]  && echo "Midnight Commander utility is not available, Install it" && exit 1

which xmgrace &> /dev/null

[ $? -ne 0 ]  && echo "Grace utility is not available, Install it" && exit 1

if [ ! -e "variables.txt" ]; then

VARTXT="liganditp = ligand.itp
ligandname = LIG
proteinpdb = protein.pdb
ligandpdb = ligand.pdb

watermodel = tip3p
forcefield = charmm27
distbox = 2
formabox = cubic

pname = NA
nname = CL"

echo "$VARTXT" >> variables.txt

fi

if [ ! -e "em.mdp" ]; then

EMMDP="integrator	= steep
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
pbc             = xyz"

echo "$EMMDP" >> em.mdp

fi

if [ ! -e "nvt.mdp" ]; then

NVTMDP="define      = -DPOSRES

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
gen-seed    = -1"

echo "$NVTMDP" >> nvt.mdp

fi

if [ ! -e "npt.mdp" ]; then

NPTMDP="define      = -DPOSRES

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
gen-vel     = no"

echo "$NPTMDP" >> npt.mdp

fi


if [ ! -e "md.mdp" ]; then

MDMDP="integrator  = md
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
gen-vel     = no"

echo "$MDMDP" >> md.mdp

fi

##############################################################################
#                      Functions definitions                                 #
##############################################################################


###### edit function: the most basic text editor: NANO #########

function edit () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68; bash"

}

###### vars function: basic configuration file #######

function vars () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 variables.txt; bash"

}

######### topols function: it will create the receptor topology ########

function topols () {

export proteinpdb=$(grep 'proteinpdb' variables.txt | sed 's/^.*= //')
export forcefield=$(grep 'forcefield' variables.txt | sed 's/^.*= //')
export watermodel=$(grep 'watermodel' variables.txt | sed 's/^.*= //')

if [ -e "$proteinpdb" ]; then
xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "gmx pdb2gmx -f $proteinpdb -ff $forcefield -water $watermodel -ignh -o protein-complex.pdb; bash"
else
dialog --title 'Atención' --msgbox 'Protein file $proteinpdb not found' 10 20

fi
}

########## edtop function #########

function edtop () {

if [ -e "topol.top" ]; then
xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 topol.top; bash"
else
dialog --title 'Atención' --msgbox 'Topology file topol.top, not found' 10 20

fi

}

########### editp function ############
function editp () {

export liganditp=$(grep 'liganditp' variables.txt | sed 's/^.*= //')

if [ -e "$liganditp" ]; then
xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 $liganditp ; bash"
else

dialog --title 'Atención' --msgbox 'El archivo de topología del ligando, no se encuentra' 10 20

fi

}

###### neutralise function: crea la caja, solvata y neutraliza ####

function neutral () {

export distbox=$(grep 'distbox' variables.txt | sed 's/^.*= //')
export formabox=$(grep 'formabox' variables.txt | sed 's/^.*= //')
export pname=$(grep 'pname' variables.txt | sed 's/^.*= //')
export nname=$(grep 'nname' variables.txt | sed 's/^.*= //')


xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e " \
gmx editconf -f protein-complex.pdb -o protein-complex-box.pdb -c -d $distbox -bt $formabox && \
gmx solvate -cs -cp protein-complex-box.pdb -o protein-complex-solv.pdb -p topol.top && \
gmx grompp -f em.mdp -c protein-complex-solv.pdb -p topol.top -o ions.tpr && \
gmx genion -s ions.tpr -o protein-complex-neutral.pdb -p topol.top -pname $pname -nname $nname -neutral \
; bash"

}

######### function para editar .mdp ##########

function edem () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 em.mdp; bash"

}

function ednvt () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 nvt.mdp; bash"

}

function ednpt () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 npt.mdp; bash"

}

function edmd () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e "nano -wcSr68 md.mdp; bash"

}

########### corridas function ############

#EM


function emrun () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e " \
gmx grompp -f em.mdp -c protein-complex-neutral.pdb -p topol.top -o em.tpr && \
gmx mdrun -v -deffnm em ; bash"

}

#NVT

function nvtrun () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e " \
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr && \
gmx mdrun -v -deffnm nvt ; bash"

}

#NPT

function nptrun () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e " \
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr && \
gmx mdrun -v -deffnm npt ; bash"

}

#MD

function mdrun () {

xterm -rightbar -bg white -fg black -fa 'Monospace' -fs 12 -geometry 92x36 -e " \
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr && \
gmx mdrun -v -deffnm md ; bash"

}



##############################################################################
#                          MENU DE LA APLICACION                             #
##############################################################################

trap 'rm menuchoices.*'  EXIT     # borrar el archivo temporal

while :
do

# Dialog utility to display options list

    dialog --clear --backtitle "Interfaz de GROMACS" --title "Menú principal" \
    --menu "Use [UP/DOWN] key to move" 20 70 35 \
    "EDIT"   "Editor básico de texto (nano)" \
    "VARS"   "Editor del archivo de variables (variables.txt)" \
    "TOPOLS" "Creador de topología de la proteína (topol.top)" \
    "EDTOP"  "Editor del archivo de topología (topol.top)" \
    "EDITP"  "Editor de la topología (.itp) del ligando" \
    "NEUTRAL" "Solvatación y neutralización" \
    "EDEM"   "Editor de los parámetros para EM" \
    "EM"     "Minimización" \
    "EDNVT"  "Editor de los parámetros para NVT" \
    "NVT"    "Equilibrado a T constante" \
    "EDNPT"  "Editor de los parámetros para NPT" \
    "NPT"    "Equilibrado a P constante" \
    "EDMD"   "Editor de los parámetros para MD" \
    "MD"     "Corrida de dinámica" \
    "EXIT"   "Salir" 2> menuchoices.$$

    retopt=$?
    choice=`cat menuchoices.$$`

    case $retopt in

           0) case $choice in

                  EDIT)   edit ;;
				  VARS)   vars ;;
                  TOPOLS) topols ;;
                  EDTOP)  edtop ;;
                  EDITP)  editp ;;
                  NEUTRAL) neutral ;;
                  EDEM)   edem ;;
                  EDNVT)  ednvt ;;
                  EDNPT)  ednpt ;;
                  EDMD)   edmd ;;
                  EM)     emrun ;;
                  NVT)    nvtrun ;;
                  NPT)    nptrun ;;
                  MD)     mdrun ;;
                  EXIT)   clear; exit 0;;

              esac ;;

          *)clear ; exit ;;
    esac

done
