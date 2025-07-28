#!/bin/bash

gcc lbm.c demo1_massdiffusion.c -o demo1 -lSDL2 -lm
gcc lbm.c demo2_cavityflow.c -o demo2 -lSDL2 -lm

#./demo1
./demo2