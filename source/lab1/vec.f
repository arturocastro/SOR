C
C MEMORY HIERARCHY DETECTION CODE
C THIS VERSION WRITTEN BY MARK BULL 8/9/97
c
	integer i,j,n,m,nops,npts,maxlength,stride,length
	parameter (nops=50000000,npts=100,maxlength=10000000)
	real x(maxlength)
        real lower,upper
	real time,junk(2)
        real res
        common /heap/ x 
c
C
C INITIALISE VECTOR DATA
C
        do i=1,maxlength
          x(i) =  3.142
        end do
        res=0.0
c
        print *, "Length     Time     Mflop/s"
c
C  USE NON-UNIT STRIDE TO REDUCE REUSE OF CACHE LINES
        stride=16
C  USE NPTS POINTS EQUALLY SPACED (LOGARITHMICALLY) BETWEEN LOWER
C  UPPER BOUNDS
        lower=100.
        upper=real(maxlength)
        do j=1,npts
           n=int(exp(log(lower)+log(upper/lower)*
     $           real(j-1)/real(npts-1)))
           m=nops/n 
C START TIMING
           time= etime(junk)
C REPEATED VECTOR SUMS
           do i=1,m*stride
              do k=1,n,stride
                 res=res+x(k)
              end do
           end do
C END TIMING 
           time=etime(junk)-time
C PREVENT DEAD CODE ELIMINATION!
           if (res .lt. 0.0) print *,res
C
C  PRINT RESULTS
           length = 4*n
           perf = real(n*m)/(1000000.*time)
           print *, length, time, perf
        end do
c
	end
