#!/bin/bash
#This file will generate the force set from phonopy_disp.yaml and {vasprun.xml}

#Directory that contain the set of calculation with different displacement; calculation file names are "disp-XXX" with XXX being the number of generated displacement;
DIR_force="../1_Run_forces_re"

#copy the phonopy_disp.yaml
cp ../0_Generate_inputs/phonopy_disp.yaml .

#copy the input: 
cp ./pp_inputs/*in . 

#activate the conda env and run phonopy
PHO=~/miniconda3/envs/west/bin/phonopy #lux

#get the series of files:
DIR=$(ls $DIR_force|grep disp)

#a strong of file list
Str_disp=''

for dir in $DIR
do
#if the file exist, append to string; otherwise end the problem. One should go back to check if the vasp calculation finished.
if [ -e "$DIR_force/$dir/scf.out" ];then
  Str_disp+="$DIR_force/$dir/scf.out "
else
  echo "no scf.out in $dir, please double check the calculation"
  exit 1
fi
done

echo found series of vasp file: $Str_disp



###If start from QE, also require scf.in ; where geometry is the one without displacement;
cp ../0_Generate_inputs/scf.in .


$PHO -f $Str_disp

