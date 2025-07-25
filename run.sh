#!/bin/bash

gcc lbm.c demo1_cavityflow.c -o cavityflow -lSDL2 -lm

./cavityflow