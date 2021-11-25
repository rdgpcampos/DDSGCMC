#!/bin/bash -x

set -eu

LAMMPS_SRC=$HOME/lammps-stable_12Dec2018/lammps-stable_12Dec2018/src

## DD-SGCMC
#MYSRC=gcmc_MD
#BIN=sgcmc

## Database construnction
MYSRC=dbcon
BIN=db

# default
mpicxx -I${LAMMPS_SRC} -c ${MYSRC}.cpp -std=c++11 -O2
mpicxx -L${LAMMPS_SRC} ${MYSRC}.o -llammps -std=c++11 -O2 -o $BIN

# intel mpi
#mpiicpc -I${LAMMPS_SRC} -c ${MYSRC}.cpp -std=c++11 -O2
#mpiicpc -L${LAMMPS_SRC} ${MYSRC}.o -llammps -std=c++11 -O2 -o $BIN
