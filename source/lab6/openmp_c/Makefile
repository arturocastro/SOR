MPIDIR = /usr/local/bin
#MPIDIR = /usr/local/software/mpich2-install/bin
CC = gcc
CCFLAGS = -O3 -Wall -Wextra -lm -fopenmp -std=c99
F90 = $(MPIDIR)/mpif90
OMP_FLAG = -fopenmp
FFLAGS = -O5 ${OMP_FLAG}
INCS = -I/local/COMP60611/utils/include/ia64/ifort
LFLAGS = -L/local/COMP60611/utils/lib/ia64/ifort -ltimers

.SUFFIXES : .o .f90 .c

sorrb: sorrb.c
	$(CC) $(CCFLAGS) -o $@${SIZE} $@.c -DM=${SIZE}

sor_jacobi: sor_jacobi.c
	$(CC) $(CCFLAGS) -o $@${SIZE} $@.c -DM=${SIZE}

clean:	
	rm -f sorrb sorrb.o sor_jacobi sor_jacobi.o
