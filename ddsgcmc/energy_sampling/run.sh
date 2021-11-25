#!/bin/bash

BIN=../../source/db

INP=inlammps
MCT=mc.ctl
RND_SEED=4238712

LMP_PATH=$HOME/lammps-stable_12Dec2018/lammps-stable_12Dec2018/src
export LD_LIBRARY_PATH=${LMP_PATH}:${LD_LIBRARY_PATH}

mpirun -np 8 $BIN $INP $MCT $RND_SEED 1>stdout 2>stderr
