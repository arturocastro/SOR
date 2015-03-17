      program sblas
      use timers
      type(TimerType) :: vecaddTime,vecsumTime,trimvTime,totalTime

      integer i,j,n,m,k,nsteps
      parameter (n=1024*60,m=512)
      parameter (nsteps=1000)
      real*8 a(n),b(n),x(n),y(m),l(m,m)
      common /heap/ a,b,x,y,l
      real*8 res
C
      call TimerCreate(vecaddTime,"vecadd")
      call TimerCreate(vecsumTime,"vecsum")
      call TimerCreate(trimvTime,"trimv")
      call TimerCreate(totalTime,"total")
C
c
c Loop over whole program to get reasonable execution time
c (simulates a time-steping code)
c
      call TimerOn(totalTime)
c
      do k=1,nsteps
c
c 'Update' vector data
c
      call vector_update(n,x,a,b)
c
c vector addition operation
c
      call TimerOn(vecaddTime)
      call vecadd (n,x,a,b)
      call TimerOff(vecaddTime)
c
c sum operation
c
      call TimerOn(vecsumTime)
      call vecsum (n,res,x)
      call TimerOff(vecsumTime)
c
c 'Update' lower triangular matrix and y vector
c
      call matrix_update (m,y,l)
c
c triangular matrix times vector operation
c
      call TimerOn(trimvTime)
      call trimv (m,y,l,a)
      call TimerOff(trimvTime)
c
c end time-step loop
	end do
c
      call TimerOff(totalTime)
c
c print out some results
c
      print *,'x(1) = ',x(1),'x(n) = ', x(n)
      print *,'Sum of x() = ',res
      print *,'y(1) = ',y(1),'y(m) = ', y(m)
C
C  OUTPUT TIMING DATA
C
      call TimerPrint(vecaddTime)
      call TimerPrint(vecsumTime)
      call TimerPrint(trimvTime)
      call TimerPrint(totalTime)
c
c
      end
c
c Subroutines
c
      subroutine vector_update (n,x,a,b)
c
c .. Scalar arguments
      integer n
c
c .. Array arguments
      real*8 x(n),a(n),b(n)
c
c .. Local scalars
      integer i
c
c .. External functions
c
      do i=1,n
         x(i) = 0.
         a(i) = 3.142
         b(i) = 3.142
      end do
c
      return
      end
c
      subroutine matrix_update (m,y,l)
c
c .. Scalar arguments
      integer m
c
c .. Array arguments
      real*8 y(m),l(m,m)
c
c .. Local scalars
      integer i,j
c
c .. External functions
c
      do i=1,m
         y(i) = 0.
         do j=1,i
            l(j,i) = 2.
         end do
      end do
c
      return
      end
c
      subroutine vecadd (n,x,a,b)
c
c .. Scalar arguments
      integer n
c
c .. Array arguments
      real*8 x(n),a(n),b(n)
c
c .. Local scalars
      integer i
c
c .. External functions
c
      do i=1,n
         x(i) = a(i) + b(i)
      end do
      return
      end
c
      subroutine vecsum (n,sum,a)
c
c .. Scalar arguments 
      integer n
      real*8 sum
c
c .. Array arguments
      real*8 a(n)
c
c .. Local scalars
      integer i
c
      sum = 0.0
      do i=1,n
        sum = sum +a(i)
      end do
      return
      end
c
      subroutine trimv (m,y,l,a)
c
c .. Scalar arguments
      integer m
c
c .. Array arguments
      real*8 y(m),l(m,m),a(m)
c
c .. Local scalars
      integer i,j
c
c .. External functions
c
      do i=1,m
        do j=1,i
          y(i) = y(i) + l(j,i)*a(j)
        end do
      end do
      return
      end

