#!/bin/bash

PL_CODE='<path to the PL code root directory>/code'

python3 -u $PL_CODE/main.py | tee pl.out

