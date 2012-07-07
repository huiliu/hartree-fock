# target: [pre]
# [TAB]  command
all: hf read2e 2e basis.o overlap

common.o: common.c common.h
	gcc common.c -c -Wall -g -lgsl -lm -o common.o

hf: hf.c common.o
	gcc hf.c common.o -Wall -g  -lgsl -lm -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e > 2e

basis.o: basis.c
	gcc basis.c -c -Wall -g -lm -o basis.o

overlap: overlap.c basis.o
	gcc overlap.c basis.o -Wall -g -lm -lgsl -o overlap

clean:
	rm -rf hf read2e *.o 2e *.mod

