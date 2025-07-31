#!/bin/bash

gcc lbm.c demo1_shockpoint.c -o demo1 -lSDL2 -lm
gcc lbm.c demo2_cavityflow.c -o demo2 -lSDL2 -lm
#gcc lbm.c demo3_channelflow.c -o demo3 -lSDL2 -lm

#./demo1
./demo2
#./demo3