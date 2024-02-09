#!/bin/bash

PL_CODE='/Users/szhang943/Programs/PL-lineshape-code_dev/code'

python3 -u $PL_CODE/main.py | tee pl.out

echo "Done:"; date
