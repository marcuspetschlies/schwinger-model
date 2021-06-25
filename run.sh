#!/bin/bash

export OMP_NUM_THREADS=4

HEAT=1.
M0=1.
BC=-1.
SEED=4
BETA=0.9

NITER=100
MEAS_EVERY=10
TAU=1.
NMD=8

LL=16
TT=32

# valgrind -v --show-reachable=yes --leak-check=full \

./schwinger -s $SEED -i $HEAT -m $M0 -c $BC -N $NITER -e $MEAS_EVERY -t $TAU -M $NMD -L $LL -T $TT 1> out 2>err

