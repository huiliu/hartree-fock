# target: [pre]
# [TAB]  command
CC = icc
OPT = -O3 -lm -lgsl -fopenmp
#OPT = -Wall -g -lm -lgsl -fopenmp -pg

all: hf read2e 2e basis.o overlap.o common.o hamiltonian.o int2e.o int scf

common.o: common.c common.h
	$(CC) common.c -c $(OPT) -o common.o

hf: hf.c common.o
	$(CC) hf.c common.o $(OPT) -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e 4 x2e.int> 2e

basis.o: basis.c basis.h
	$(CC) basis.c -c $(OPT) -o basis.o

overlap.o: overlap.c overlap.h
	$(CC) overlap.c $(OPT) -c -o overlap.o

hamiltonian.o: hamiltonian.c hamiltonian.h
	$(CC) hamiltonian.c $(OPT) -c -o hamiltonian.o

int2e.o: int2e.c int2e.h
	$(CC) int2e.c -c $(OPT) -o int2e.o

int: int.c common.o basis.o overlap.o hamiltonian.o int2e.o 
	$(CC) int.c common.o basis.o overlap.o hamiltonian.o int2e.o $(OPT) -o int

scf: scf.c scf.h common.o basis.o overlap.o hamiltonian.o int2e.o
	$(CC) scf.c basis.o overlap.o common.o hamiltonian.o int2e.o $(OPT) -o scf

clean:
	rm -rf hf read2e *.o 2e *.mod int scf
