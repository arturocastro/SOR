program main
  use timers
  type(TimerType) :: gbTime
  integer check
  
  integer maxn,n,i
  parameter (n=20000)
  integer numpairs(n)

  call TimerCreate(gbTime,"goldbach")

  do i=1,n
     numpairs(i) = 0 
  end do

! FIND NO. OF GOLDBACH PAIRS FOR THE FIRST n EVEN NUMBERS

  call TimerOn(gbTime)	
  do i=1,n
     call goldbach(2*i,numpairs(i))
  end do
  call TimerOff(gbTime)

!  OUTPUT SELECTED RESULTS

  do i=1,n,n/10
     print *, 2*i, numpairs(i)
  end do

!  OUTPUT TIMING DATA

  call TimerPrint(gbTime)

  stop 
end program main

subroutine goldbach(evenno,numpairs)

  integer evenno,numpairs,j
  logical isprime


! TEST ALL PAIRS WHICH SUM TO evenno FOR PRIMALITY
	
  do j=2,evenno/2
     if (isprime(j).and.isprime(evenno-j)) numpairs = numpairs+1
  end do

  return 
end subroutine goldbach

logical function isprime(num)

  integer num,i


! CRUDE PRIMALITY TEST: CHECK ALL INTEGERS LESS THAN SQARE ROOT
! OF NUMBER 
 	
  isprime = .true.
  i=2 
  do while ( (i.le.(int(sqrt(real(num)))) ) .and. isprime)
     isprime = mod(num,i).ne.0
     i=i+1
  end do

  return
end function isprime
