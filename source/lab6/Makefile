MPIDIR = /usr/local/bin
#MPIDIR = /usr/local/software/mpich2-install/bin
F90 = $(MPIDIR)/mpif90
OMP_FLAG = -fopenmp
FFLAGS = -O5 ${OMP_FLAG}
INCS = -I/local/COMP60611/utils/include/ia64/ifort
LFLAGS = -L/local/COMP60611/utils/lib/ia64/ifort -ltimers

.SUFFIXES : .o .f90

sorf: sorf.f90
	$(F90) $(FFLAGS) -o $@ $@.f90 ${INCS} ${LFLAGS}

sorrbf: sorrbf.f90
	$(F90) $(FFLAGS) -o $@ $@.f90 ${INCS} ${LFLAGS}
clean:	
	rm -f sorf
