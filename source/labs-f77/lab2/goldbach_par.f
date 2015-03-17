      program main
      use timers
      type(TimerType) :: gbTime

      integer maxn,n,i
      parameter (n=20000)
      integer numpairs(n)
C
      call TimerCreate(gbTime,"goldbach")
C
      do i=1,n
         numpairs(i) = 0 
      end do      
C
C FIND NO. OF GOLDBACH PAIRS FOR THE FIRST n EVEN NUMBERS
C
      call TimerOn(gbTime)
!$omp parallel do default (none) private (i)
!$omp+  shared (numpairs) schedule(static)
      do i=1,n
         call goldbach(2*i,numpairs(i))
      end do
      call TimerOff(gbTime)
C
C  OUTPUT SELECTED RESULTS
C
      do i=1,n,n/10
         print *, 2*i, numpairs(i)
      end do
C
C  OUTPUT TIMING DATA
C
      call TimerPrint(gbTime)
C
      stop 
      end
C
      subroutine goldbach(evenno,numpairs)
c
      integer evenno,numpairs,j
      logical isprime
c
C
C TEST ALL PAIRS WHICH SUM TO evenno FOR PRIMALITY
C	
      do j=2,evenno/2
         if (isprime(j).and.isprime(evenno-j)) numpairs = numpairs+1
      end do
c
      return 
      end 
c
c
c
      logical function isprime(num)
c
      integer num,i
c
C
C CRUDE PRIMALITY TEST: CHECK ALL INTEGERS LESS THAN SQARE ROOT
C OF NUMBER 
C 	
      isprime = .true.
      i=2 
      do while ( (i.le.(int(sqrt(real(num)))) ) .and. isprime)
         isprime = mod(num,i).ne.0
         i=i+1
      end do
c
      return
      end
