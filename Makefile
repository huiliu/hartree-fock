# target: [pre]
# [TAB]  command
CC = gcc
#OPT = -O3 -lm -lgsl -fopenmp
OPT = -Wall -g -lm -lgsl -fopenmp

all: hf read2e 2e basis.o overlap.o common.o hamiltonian.o ints.o int2e.o int_s test scf

common.o: common.c common.h
	$(CC) common.c -c $(OPT) -o common.o

hf: hf.c common.o
	$(CC) hf.c common.o $(OPT) -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e 4 x2e.int> 2e

ints.o: ints.c
	$(CC) ints.c   $(OPT) -c -o ints.o

basis.o: basis.c basis.h
	$(CC) basis.c -c $(OPT) -o basis.o

overlap.o: overlap.c overlap.h
	$(CC) overlap.c $(OPT) -c -o overlap.o

coord: coord.c basis.o overlap.o common.o
	$(CC) coord.c basis.o overlap.o common.o $(OPT) -o coord

hamiltonian.o: hamiltonian.c hamiltonian.h
	$(CC) hamiltonian.c $(OPT) -c -o hamiltonian.o

int2e.o: int2e.c int2e.h
	$(CC) int2e.c -c $(OPT) -o int2e.o

int_s: int_s.c basis.o common.o overlap.o ints.o
	$(CC) int_s.c basis.o overlap.o common.o ints.o $(OPT) -o int_s

test: test.c common.o basis.o overlap.o hamiltonian.o int2e.o ints.o
	$(CC) test.c common.o basis.o overlap.o hamiltonian.o int2e.o ints.o $(OPT) -o test

scf: scf.c scf.h common.o basis.o overlap.o hamiltonian.o int2e.o
	$(CC) scf.c basis.o overlap.o common.o hamiltonian.o ints.o int2e.o $(OPT) -o scf

clean:
	rm -rf hf read2e *.o 2e *.mod hamiltonian coord ints test int_s scf
