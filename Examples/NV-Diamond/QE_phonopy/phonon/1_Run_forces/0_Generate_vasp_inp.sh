#!/bin/bash

#Directory of poscar
Dir_POS="../0_Generate_inputs/"

#Directory include qe inputs header
QE_Sample="./Sample_input"

POS_file=`ls $Dir_POS | grep supercell-`

for f_pos in ${POS_file[@]}
do
#grep the number of displacement card:
number=$(echo $f_pos | grep -o '[0-9]*')


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


#copy input from sample file to displacement file:
cp $QE_Sample/* $filename

#copy supercell input to the displacement file
cp $Dir_POS/$f_pos ./$filename/

#append the coordinate to scf.in header file in Sample_input:
cd $filename
cat header.in $f_pos >| scf.in
cd ..
done
