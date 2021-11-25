#!/bin/bash

BIN=../../source/dd

INP=inlammps
N=1
MCT=mc.ctl
RND_SEED=5476378

LMP_PATH=$HOME/lammps-stable_12Dec2018/lammps-stable_12Dec2018/src
export LD_LIBRARY_PATH=${LMP_PATH}:${LD_LIBRARY_PATH}

mpirun -np 8  $BIN $N $INP $MCT 1>stdout 2>stderr

