        program main
        use timers
        type(TimerType) :: seriesTime,monteTime,totalTime

        integer npoints,nthreads,maxthreads,i,ierr
        parameter (npoints=300000,maxthreads=16)
        complex*16 c(npoints)
        real*8 r1,r2,halfpi,err1,area,err2
        real rand
        integer numinside,lnuminside(maxthreads)
        integer omp_get_num_threads
        external monte
        external omp_set_num_threads,omp_get_num_threads
C
        call TimerCreate(seriesTime,"series")
        call TimerCreate(monteTime,"monte")
        call TimerCreate(totalTime,"Series+monte") 
C
C GENERATE RANDOM NUMBERS
C
      call srand(54321)
      do i=1,npoints
       	 r1=rand()
         r2=rand()
         c(i)=CMPLX(-2.0+2.5*r1,1.125*r2)
      end do
C
      call TimerOn(totalTime)
C
C  CALCULATE PI/2 FROM GREGORY'S FORMULA
C
      call TimerOn(seriesTime)
      call series(halfpi,err1)
      call TimerOff(seriesTime)
C
C  CALCULATE AREA OF MANDELBROT SET BY MONTE CARLO SAMPLING
C
      call TimerOn(monteTime)
C  INVOKE monte in parallel
!$omp parallel
C  GET NUMBER OF THREADS BEING USED
      nthreads = omp_get_num_threads()
C
      call monte(c,lnuminside)
!$omp end parallel
C
C  SUM UP CONTRIBUTIONS FROM EACH THREAD 
C	
      numinside = 0
      do i=1,nthreads
         numinside = numinside + lnuminside(i)
      end do
C
      call TimerOff(monteTime)
      call TimerOff(totalTime)
C
C  OUTPUT RESULTS
C
      area = 2.0*2.5*1.125 * real(numinside)/real(npoints)
      err2 = area/sqrt(real(npoints))
      print *, "Pi/2 = ", halfpi," +/- ",err1
      print *, "Area of Mandelbrot set = ",area," +/- ",err2
C
C
C  OUTPUT TIMING DATA
C
        call TimerPrint(seriesTime)
        call TimerPrint(monteTime)
        call TimerPrint(totalTime)
C
      stop
      end
C
      subroutine series(halfpi,err)
c
      real*8 halfpi,err,sum
      integer sign,terms,denom,i
      parameter (terms = 100000000)
c
      sum=0.0d0
      sign=1
      denom=1
c
      do i=1,terms
         sum=sum + sign*(1.0d0/denom)
         sign=(-sign)
         denom=denom + 2
      end do
c
      halfpi = 2.0d0 * sum 
      err = 2.0d0/terms
      return
      end 
	
      subroutine monte(c,lnuminside)
c
      integer npoints,maxiter,i,j,num,slice,myid,ilo,ihi
      integer nthreads,maxthreads
      parameter (npoints=300000,maxiter=10000,maxthreads=16)
      complex*16 c(npoints),z
      integer lnuminside(maxthreads)
      integer omp_get_num_threads,omp_get_thread_num
      external omp_get_num_threads,omp_get_thread_num
c
      nthreads = omp_get_num_threads()
      slice = (npoints+nthreads-1)/nthreads
      myid =  omp_get_thread_num()
      ilo = slice * myid + 1
      ihi = min(npoints,(myid+1)*slice)

      num=0
      do i=ilo,ihi
         z=(0.0,0.0)
         do j=1,maxiter
            z = z*z +c(i)	 
            if (abs(z).gt.2.0) goto 10
         end do
         num = num +1
 10      continue
      end do 
C
      lnuminside(myid+1) = num
C
      return
      end	



