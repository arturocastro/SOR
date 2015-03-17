FC = gfortran
OMP_FLAG = -fopenmp
FFLAGS = -O -c ${OMP_FLAG}
INCS = -I/localhome/COMP60611/utils/include/ia64/gfortran
LFLAGS = -L/localhome/COMP60611/utils/lib/ia64/gfortran -ltimers

.f.o:
	${FC} ${FFLAGS} -o $@ $*.f ${INCS}

sblas:sblas.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

sblas_par:sblas_par.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

sblas_par2:sblas_par2.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o sblas sblas_par sblas_par2 *.o

