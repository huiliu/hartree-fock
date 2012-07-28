# target: [pre]
# [TAB]  command
all: hf read2e 2e basis.o overlap.o common.o hamiltonian.o ints.o int2e.o check_2e int_s test scf

common.o: common.c common.h
	gcc common.c -c -Wall -g -lgsl -lm -o common.o

hf: hf.c common.o
	gcc hf.c common.o -Wall -g  -lgsl -lm -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e > 2e

ints.o: ints.c
	gcc ints.c -lm -Wall -g -c -o ints.o

basis.o: basis.c basis.h
	gcc basis.c -c -Wall -g -lm -lgsl -o basis.o

overlap.o: overlap.c overlap.h basis.o common.o
	gcc overlap.c basis.o common.o -c -Wall -g -lm -lgsl -o overlap.o

coord: coord.c basis.o overlap.o common.o
	gcc coord.c basis.o overlap.o common.o -Wall -g -lgsl -lm -o coord

hamiltonian.o: hamiltonian.c hamiltonian.h basis.o overlap.o common.o ints.o
	gcc hamiltonian.c basis.o overlap.o common.o ints.o -c -Wall -g -lgsl -lm -o hamiltonian.o

int2e.o: int2e.c int2e.h basis.o common.o overlap.o
	gcc int2e.c basis.o overlap.o common.o -c -Wall -g -lgsl -lm -o int2e.o

check_2e: check_2e.c basis.o common.o overlap.o ints.o
	gcc check_2e.c basis.o overlap.o common.o ints.o -Wall -g -lgsl -lm -o check_2e

int_s: int_s.c basis.o common.o overlap.o ints.o
	gcc int_s.c basis.o overlap.o common.o ints.o -Wall -g -lgsl -lm -o int_s

test: test.c basis.o common.o overlap.o hamiltonian.o ints.o int2e.o
	gcc test.c basis.o overlap.o common.o hamiltonian.o ints.o int2e.o -Wall -g -lgsl -lm -o test

scf: scf.c scf.h common.o basis.o overlap.o hamiltonian.o int2e.o
	gcc scf.c basis.o overlap.o common.o hamiltonian.o ints.o int2e.o -Wall -g -lgsl -lm -o scf

clean:
	rm -rf hf read2e *.o 2e *.mod hamiltonian coord ints test int_s check_2e
