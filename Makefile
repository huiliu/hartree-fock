# target: [pre]
# [TAB]  command
CC = icc
OPT = -O3 -lm -lgsl -fopenmp -pg
#OPT = -Wall -g -lm -lgsl -fopenmp -pg

all: basis.o overlap.o common.o hamiltonian.o int2e.o int1e.o eri.o int scf redblack.o

common.o: common.c common.h
	$(CC) common.c -c $(OPT) -o common.o

hf: hf.c common.o redblack.o
	$(CC) hf.c common.o redblack.o $(OPT) -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e 4 x2e.int> 2e

basis.o: basis.c basis.h
	$(CC) basis.c -c $(OPT) -o basis.o

redblack.o: redblack.c redblack.h
	$(CC) redblack.c -c $(OPT) -o redblack.o

overlap.o: overlap.c overlap.h
	$(CC) overlap.c $(OPT) -c -o overlap.o

hamiltonian.o: hamiltonian.c hamiltonian.h
	$(CC) hamiltonian.c $(OPT) -c -o hamiltonian.o

int1e.o: int1e.c int1e.h
	$(CC) int1e.c -c $(OPT) -o int1e.o

eri.o: eri.c eri.h
	$(CC) eri.c -c $(OPT) -o eri.o

int2e.o: int2e.c int2e.h
	$(CC) int2e.c -c $(OPT) -o int2e.o

int: int.c common.o basis.o overlap.o hamiltonian.o int1e.o int2e.o eri.o redblack.o
	$(CC) int.c common.o basis.o overlap.o hamiltonian.o int1e.o int2e.o eri.o redblack.o $(OPT) -o int

scf: scf.c scf.h common.o basis.o overlap.o hamiltonian.o int1e.o int2e.o eri.o redblack.o
	$(CC) scf.c basis.o overlap.o common.o hamiltonian.o int2e.o int1e.o eri.o redblack.o $(OPT) -o scf

clean:
	rm -rf hf read2e *.o 2e *.mod int scf gmon.out
