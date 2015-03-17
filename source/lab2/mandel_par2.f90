program main
  use timers
  type(TimerType) :: seriesTime,monteTime,totalTime
  
  integer npoints,i,numinside
  parameter (npoints=300000)
  complex*16 c(npoints)
  real*8 r1,r2,halfpi,err1,area,err2
  real rand

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
  call monte(c,numinside)
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
	
subroutine monte(c,numinside)

  integer npoints,npoints_p,maxiter,maxiter_p,numinside,i,j
  parameter (npoints_p=300000,maxiter_p=10000)
  complex*16 c(npoints_p),z

  maxiter=maxiter_p
  npoints=npoints_p
  numinside=0
!$omp parallel do default (none) private(i,j,z) &
!$omp shared (c,maxiter,npoints) reduction(+:numinside)

  do i=1,npoints
     z=(0.0,0.0)
     do j=1,maxiter
        z = z*z +c(i)	 
        if (abs(z).gt.2.0) goto 10
     end do
     numinside = numinside +1
10   continue
  end do

  
  return
end subroutine monte

