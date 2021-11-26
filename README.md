# DDSGCMC
C++ source codes and test cases for the database-driven semi-grand canonical Monte Carlo simulation method.

If you find this work useful in your research, please cite:

[R.P. Campos, S. Shinzato, A. Ishii, S. Nakamura, S. Ogata, Database-driven semi-grand canonical Monte Carlo method:
application to segregation isotherm on defects in alloys, Physical Review E **104**, 025310 (2021)](https://doi.org/10.1103/PhysRevE.104.025310)


The source codes used in this method can be accessed in the source directory, while the run scripts and all necessary files are
located in the ddsgcmc directory. Note that to run the NiCoCr example, first the NiCoCr.lammps.eam potential
needs to be included in the ./ddsgcmc directory (See the README file in the same directory).

These scripts were tested with the stable 12th of December of 2018 version of LAMMPS. Before using them, modify the path to
the LAMMPS src files in source/Makefile (INC_PATH and LIB_PATH). It is also necessary to modify the same path in the run.sh scripts
in /ddsgcmc/dd_sgcmc and /ddsgcmc/energy_sampling before running either of these calculations.

Contact the authors for more information or comments regarding this method:

rdgpcampos@tsme.me.es.osaka-u.ac.jp
