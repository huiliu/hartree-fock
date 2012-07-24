# target: [pre]
# [TAB]  command
all: hf read2e 2e basis.o overlap.o common.o kinetic.o nuclear_elect.o ints.o int2e.o check_2e int_s test

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

kinetic.o: kinetic.c kinetic.h overlap.o basis.o common.o ints.o
	gcc kinetic.c overlap.o basis.o common.o ints.o -c -Wall -g -lgsl -lm -o kinetic.o

coord: coord.c basis.o overlap.o common.o
	gcc coord.c basis.o overlap.o common.o -Wall -g -lgsl -lm -o coord

nuclear_elect.o: nuclear_elect.c nuclear_elect.h basis.o overlap.o common.o ints.o
	gcc nuclear_elect.c basis.o overlap.o common.o ints.o -c -Wall -g -lgsl -lm -o nuclear_elect.o

int2e.o: int2e.c int2e.h basis.o common.o overlap.o
	gcc int2e.c basis.o overlap.o common.o -c -Wall -g -lgsl -lm -o int2e.o

check_2e: check_2e.c basis.o common.o overlap.o ints.o
	gcc check_2e.c basis.o overlap.o common.o ints.o -Wall -g -lgsl -lm -o check_2e

int_s: int_s.c basis.o common.o overlap.o ints.o
	gcc int_s.c basis.o overlap.o common.o ints.o -Wall -g -lgsl -lm -o int_s

test: test.c basis.o common.o overlap.o kinetic.o nuclear_elect.o ints.o int2e.o
	gcc test.c basis.o overlap.o common.o kinetic.o nuclear_elect.o ints.o int2e.o -Wall -g -lgsl -lm -o test

clean:
	rm -rf hf read2e *.o 2e *.mod kinetic coord ints nuclear_elect test int_s check_2e
