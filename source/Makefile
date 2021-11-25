# DB sampling

# Need to adapt lammps src file path



PROGRAM_DB = dbcon
SRC_DB = db

PROGRAM_DD = ddsgcmc
SRC_DD = dd


MPICXX = mpicxx
INC_PATH = $(HOME)/lammps-stable_12Dec2018/lammps-stable_12Dec2018/src
LIB_PATH = $(HOME)/lammps-stable_12Dec2018/lammps-stable_12Dec2018/src
CFLAGS  = -O2 -std=c++11
LDFLAGS = -llammps

.SUFFIXES: .cpp .o

all: db dd

$(SRC_DB): $(PROGRAM_DB)
	$(MPICXX) -L$(LIB_PATH) $^.o $(LDFLAGS) $(CFLAGS) -o $@

$(PROGRAM_DB): 
	$(MPICXX) -I$(INC_PATH) -c $@.cpp $(CFLAGS)

$(SRC_DD): $(PROGRAM_DD)
	$(MPICXX) -L$(LIB_PATH) $^.o $(LDFLAGS) $(CFLAGS) -o $@

$(PROGRAM_DD): 
	$(MPICXX) -I$(INC_PATH) -c $@.cpp $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(PROGRAM) *.o
