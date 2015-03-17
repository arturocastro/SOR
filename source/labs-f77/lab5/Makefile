MPIDIR = /usr/local/software/mpich2-install/bin
F77 = $(MPIDIR)/mpif90
OMP_FLAG = -fopenmp
FFLAGS = -O0 ${OMP_FLAG}
INCS = -I/localhome/COMP60611/utils/include/ia64/gfortran
LFLAGS = -L/localhome/COMP60611/utils/lib/ia64/gfortran -ltimers

sblas_1: sblas_1.f
	$(F77) $(FFLAGS) -o $@ $@.f ${INCS} ${LFLAGS}

sblas_2: sblas_2.f
	$(F77) $(FFLAGS) -o $@ $@.f ${INCS} ${LFLAGS}

clean:	
	rm -f sblas_1 sblas_2
