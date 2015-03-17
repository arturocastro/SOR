FC = gfortran
OMP_FLAG = -fopenmp
FFLAGS = -O -c ${OMP_FLAG}
INCS = -I/localhome/COMP60611/utils/include/ia64/gfortran
LFLAGS = -L/localhome/COMP60611/utils/lib/ia64/gfortran -ltimers

.f.o:
	${FC} ${FFLAGS} -o $@ $*.f ${INCS}

aquad: aquad.o
	${FC} ${OMP_FLAG} -o $@ $@.o ${LFLAGS}

clean:
	rm *.o aquad aquad_par 
