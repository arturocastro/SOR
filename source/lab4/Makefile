FC = gfortran
OMP_FLAG = -fopenmp
FFLAGS = -O -c ${OMP_FLAG}
INCS = -I/local/COMP60611/utils/include/ia64/gfortran
LFLAGS = -L/local/COMP60611/utils/lib/ia64/gfortran -ltimers

.SUFFIXES : .o .f90

.f90.o:
	${FC} ${FFLAGS} -o $@ $*.f90 ${INCS}

aquad: aquad.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o aquad aquad_par 
