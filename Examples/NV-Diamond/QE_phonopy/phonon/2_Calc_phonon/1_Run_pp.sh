#!/bin/bash

echo "before run":
ls
PHO=~/miniconda3/envs/west/bin/phonopy

#with mesh.conf, plot density of state (DOS); This only works for q mesh >1
#$PHO -p mesh.conf

#Thermo properties:
#$PHO -t mesh.conf
#plot by 
#$PHO -t -p mesh.conf

#Projected DOS: require pdos.conf input
#$PHO -p pdos.conf

#$PHO phonopy.in -v

#Run phonon band from vasp
#$PHO --qpoints="0 0 0" --writedm phonopy.in -v  

#Run phonon band from qe:
#$PHO --qe -c scf.in --qpoints="0 0 0" --writedm phonopy.in -v 
$PHO --qe -c scf.in phonopy.in -v
#Print dynamic matrix
$PHO --qe -c scf.in --qpoints="0 0 0" --writedm phonopy.in -v

echo "phonopy done"
ls
