MPIDIR = /usr/local/bin
#MPIDIR = /usr/local/software/mpich2-install/bin
CC = mpicc
CCFLAGS = -O3 -Wall -Wextra -lm -std=c99
INCS = -I/local/COMP60611/utils/include/ia64/ifort
LFLAGS = -L/local/COMP60611/utils/lib/ia64/ifort -ltimers

.SUFFIXES : .o .f90 .c

sor_jacobi: sor_jacobi.c
	$(CC) $(CCFLAGS) -o $@ $@.c

sorrb: sorrb.c
	$(CC) $(CCFLAGS) -o $@ $@.c
clean:	
	rm -f sorrb sorrb.o sor_jacobi sor_jacobi.o
