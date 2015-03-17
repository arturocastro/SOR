program sblas
  use timers
  type(TimerType) :: vecaddTime,vecsumTime,trimvTime,totalTime
  
  integer i,j,n,m,k,nsteps
  parameter (n=1024*60,m=512)
  parameter (nsteps=1000)
  real*8 a(n),b(n),x(n),y(m),l(m,m)
  common /heap/ a,b,x,y,l
  real*8 res

  call TimerCreate(vecaddTime,"vecadd")
  call TimerCreate(vecsumTime,"vecsum")
  call TimerCreate(trimvTime,"trimv")
  call TimerCreate(totalTime,"total")

! Loop over whole program to get reasonable execution time
! (simulates a time-steping code)

  call TimerOn(totalTime)

  do k=1,nsteps

! 'Update' vector data

    call vector_update(n,x,a,b)

! vector addition operation

    call TimerOn(vecaddTime)
    call vecadd (n,x,a,b)
    call TimerOff(vecaddTime)

! sum operation

    call TimerOn(vecsumTime)
    call vecsum (n,res,x)
    call TimerOff(vecsumTime)

! 'Update' lower triangular matrix and y vector

    call matrix_update (m,y,l)

! triangular matrix times vector operation

    call TimerOn(trimvTime)
    call trimv (m,y,l,a)
    call TimerOff(trimvTime)

! end time-step loop
 end do

 call TimerOff(totalTime)

! print out some results

 print *,'x(1) = ',x(1),'x(n) = ', x(n)
 print *,'Sum of x() = ',res
 print *,'y(1) = ',y(1),'y(m) = ', y(m)

!  OUTPUT TIMING DATA

 call TimerPrint(vecaddTime)
 call TimerPrint(vecsumTime)
 call TimerPrint(trimvTime)
 call TimerPrint(totalTime)
 
end program sblas

! Subroutines

subroutine vector_update (n,x,a,b)

! .. Scalar arguments
  integer n

! .. Array arguments
  real*8 x(n),a(n),b(n)

! .. Local scalars
  integer i

! .. External functions

  do i=1,n
     x(i) = 0.
     a(i) = 3.142
     b(i) = 3.142
  end do

  return
end subroutine vector_update

subroutine matrix_update (m,y,l)

! .. Scalar arguments
  integer m

! .. Array arguments
  real*8 y(m),l(m,m)

! .. Local scalars
  integer i,j

! .. External functions

  do i=1,m
     y(i) = 0.
     do j=1,i
        l(j,i) = 2.
     end do
  end do

  return
end subroutine matrix_update

subroutine vecadd (n,x,a,b)
!
! .. Scalar arguments
  integer n

! .. Array arguments
  real*8 x(n),a(n),b(n)

! .. Local scalars
  integer i

! .. External functions

  do i=1,n
     x(i) = a(i) + b(i)
  end do
  return
end subroutine vecadd

subroutine vecsum (n,sum,a)

! .. Scalar arguments 
  integer n
  real*8 sum

! .. Array arguments
  real*8 a(n)

! .. Local scalars
  integer i

  sum = 0.0
  do i=1,n
     sum = sum +a(i)
  end do
  return
end subroutine vecsum

subroutine trimv (m,y,l,a)

! .. Scalar arguments
  integer m

! .. Array arguments
  real*8 y(m),l(m,m),a(m)

! .. Local scalars
  integer i,j

! .. External functions

  do i=1,m
     do j=1,i
        y(i) = y(i) + l(j,i)*a(j)
     end do
  end do
  return
end subroutine trimv

