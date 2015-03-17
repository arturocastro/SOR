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

  call TimerCreate(seriesTime,"series")
  call TimerCreate(monteTime,"monte")
  call TimerCreate(totalTime,"Series+monte") 

! GENERATE RANDOM NUMBERS

  call srand(54321)
  do i=1,npoints
     r1=rand()
     r2=rand()
     c(i)=CMPLX(-2.0+2.5*r1,1.125*r2)
  end do

  call TimerOn(totalTime)

!  CALCULATE PI/2 FROM GREGORY'S FORMULA

  call TimerOn(seriesTime)
  call series(halfpi,err1)
  call TimerOff(seriesTime)

!  CALCULATE AREA OF MANDELBROT SET BY MONTE CARLO SAMPLING

  call TimerOn(monteTime)
!  INVOKE monte in parallel

!$omp parallel

!  GET NUMBER OF THREADS BEING USED
  nthreads = omp_get_num_threads()

  call monte(c,lnuminside)
!$omp end parallel

!  SUM UP CONTRIBUTIONS FROM EACH THREAD 
  numinside = 0
  do i=1,nthreads
     numinside = numinside + lnuminside(i)
  end do

  call TimerOff(monteTime)
  call TimerOff(totalTime)

!  OUTPUT RESULTS

  area = 2.0*2.5*1.125 * real(numinside)/real(npoints)
  err2 = area/sqrt(real(npoints))
  print *, "Pi/2 = ", halfpi," +/- ",err1
  print *, "Area of Mandelbrot set = ",area," +/- ",err2

!  OUTPUT TIMING DATA

  call TimerPrint(seriesTime)
  call TimerPrint(monteTime)
  call TimerPrint(totalTime)

  stop
end program main

subroutine series(halfpi,err)
  
  real*8 halfpi,err,sum
  integer sign,terms,denom,i
  parameter (terms = 100000000)

  sum=0.0d0
  sign=1
  denom=1

  do i=1,terms
     sum=sum + sign*(1.0d0/denom)
     sign=(-sign)
     denom=denom + 2
  end do

  halfpi = 2.0d0 * sum 
  err = 2.0d0/terms
  return
end subroutine series
	
subroutine monte(c,lnuminside)

  integer npoints,maxiter,i,j,num,slice,myid,ilo,ihi
  integer nthreads,maxthreads
  parameter (npoints=300000,maxiter=10000,maxthreads=16)
  complex*16 c(npoints),z
  integer lnuminside(maxthreads)
  integer omp_get_num_threads,omp_get_thread_num
  external omp_get_num_threads,omp_get_thread_num

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
10   continue
  end do

  lnuminside(myid+1) = num

  return
end subroutine monte



