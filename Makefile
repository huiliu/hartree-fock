# target: [pre]
# [TAB]  command
all: hf read2e 2e

hf: hf.c 
	gcc hf.c -Wall -g  -lgsl -lm -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e > 2e

clean:
	rm -rf hf read2e *.o 2e

