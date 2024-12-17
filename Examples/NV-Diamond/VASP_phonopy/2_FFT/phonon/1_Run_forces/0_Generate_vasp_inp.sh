#!/bin/bash

#Directory of poscar
Dir_POS="../0_Generate_POSCAR/"

#Directory include VASP inputs
VASP_Sample="./Sample_input"

POS_file=`ls ../0_Generate_POSCAR/ | grep POSCAR-`

for f_pos in ${POS_file[@]}
do
#grep the number of displacement card:
number=$(echo $f_pos | grep -o '[0-9]*$')


#check if directory exist, and create file:
filename="disp-$number"
if [ -e "$filename" ]; then
    # File exists, remove it and create a new one
    echo "File exists, removing and creating a new one."
    rm -rf "$filename"
    mkdir "$filename"
else
    # File does not exist, create it
    echo "Creating file $filename."
    mkdir "$filename"
fi

#copy POSCAR to the displacement 
cp $Dir_POS/$f_pos ./$filename/POSCAR

#copy input from sample file to displacement file:
cp $VASP_Sample/* $filename
done
