# target: [pre]
# [TAB]  command
CC = icc
#OPT = -O3 -lm -lgsl -fopenmp
OPT = -O3 -lm -lgsl -fopenmp -pg
#OPT = -Wall -g -lm -lgsl -fopenmp -pg

all: print.o basis.o common.o int eri_os.o eri_drive.o

common.o: common.c common.h
	$(CC) common.c -c $(OPT) -o common.o

hf: hf.c common.o
	$(CC) hf.c common.o $(OPT) -o hf

read2e: read2e.f90
	ifort read2e.f90 -o read2e

2e: read2e
	./read2e 4 x2e.int> 2e

print.o: print.c print.h
	$(CC) print.c -c $(OPT) -o print.o

basis.o: basis.c basis.h
	$(CC) basis.c -c $(OPT) -o basis.o

eri_drive.o: eri_drive.c eri_drive.h
	$(CC) eri_drive.c -c $(OPT) -o eri_drive.o

eri_os.o: eri_os.c eri_os.h
	$(CC) eri_os.c -c $(OPT) -o eri_os.o

int: int.c common.o basis.o eri_os.o eri_drive.o print.o
	$(CC) int.c common.o basis.o eri_os.o eri_drive.o print.o $(OPT) -o int

scf: scf.c scf.h common.o basis.o eri_os.o print.o
	$(CC) scf.c basis.o common.o eri_os.o eri_drive.o print.o $(OPT) -o scf

clean:
	rm -rf hf read2e *.o 2e *.mod int scf gmon.out
