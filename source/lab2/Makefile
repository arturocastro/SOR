FC = gfortran
OMP_FLAG = -fopenmp
FFLAGS = -O -c ${OMP_FLAG}
INCS = -I/local/COMP60611/utils/include/ia64/gfortran
LFLAGS = -L/local/COMP60611/utils/lib/ia64/gfortran -ltimers

.SUFFIXES : .o .f90

.f90.o:
	${FC} ${FFLAGS} -o $@ $*.f90 ${INCS}

mandel:mandel.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

mandel_par:mandel_par.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

mandel_par2:mandel_par2.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

goldbach:goldbach.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

goldbach_par:goldbach_par.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o mandel mandel_par mandel_par2 goldbach goldbach_par

