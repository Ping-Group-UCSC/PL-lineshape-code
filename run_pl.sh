#!/bin/bash
ENV=/Users/szhang943/Python_Env/west/bin
PL_CODE='/Users/szhang943/Programs/PL-lineshape-code_dev/code'

$ENV/python3 -u $PL_CODE/main.py | tee pl.out

echo "Done:"; date
