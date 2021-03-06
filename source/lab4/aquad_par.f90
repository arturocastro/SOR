PROGRAM TestQuad

  use timers
  type(TimerType) :: func1Time,func2Time,func3Time
  
  REAL*8 left, right, result, eps, Quad
  REAL*8     Func1, Func2, Func3
  EXTERNAL   Func1, Func2, Func3

!     Initialise timers
  call TimerCreate(func1Time,"Func1")
  call TimerCreate(func2Time,"Func2")
  call TimerCreate(func3Time,"Func3")

  call TimerOn(func1Time)
  left  = 0.0
  right = 1.0
  eps   = 1.0E-7
  result = Quad (Func1, left, right, eps)
  call TimerOff(func1Time)

  PRINT *, "Result for Func1 = ", result
  call TimerPrint(func1Time)

  call TimerOn(func2Time)
  left  = 0.0
  right = 1.0
  eps   = 1.0E-7
  result = Quad (Func2, left, right, eps)
  call TimerOff(func2Time)

  PRINT *, "Result for Func2 = ", result
  call TimerPrint(func2Time)

  call TimerOn(func3Time)
  left  = 0.0
  right = 1.0
  eps   = 1.0E-7
  result = Quad (Func3, left, right, eps)
  call TimerOff(func3Time)

  PRINT *, "Result for Func3 = ", result
  call TimerPrint(func3Time)

END PROGRAM TestQuad

REAL*8 FUNCTION Quad (Func, left, right, epsilon)
  REAL*8 left, right, epsilon, Func

! *** Adaptive Quadrture

  INTEGER STACKSIZE,IL,IR,IE
  PARAMETER (STACKSIZE = 1000,IL=1,IR=2,IE=3)

  REAL*8 stack  (3,STACKSIZE)
  INTEGER stackpointer

  REAL*8 result, abserror, l, r, m, est1, est2, eps

  stackpointer = 1
  stack(il,stackpointer) = left
  stack(ir,stackpointer) = right
  stack(ie,stackpointer) = epsilon

  result = 0.0

  call Stackops(stackpointer,stack,result,Func) 

  Quad = result

  RETURN
END FUNCTION Quad

SUBROUTINE Stackops(stackpointer,stack,result,Func)

  INTEGER STACKSIZE,IL,IR,IE
  PARAMETER (STACKSIZE = 1000,IL=1,IR=2,IE=3)

  REAL*8 stack(3,STACKSIZE)
  INTEGER stackpointer

  REAL*8 result, Func
  REAL*8 abserror, l, r, m, est1, est2, eps

! *** ours
  INTEGER active_threads
  active_threads = 0
! *** ours

! *** ours
!$OMP parallel default(none) shared(result, stack, stackpointer, active_threads) private(abserror, l, r, m, est1, est2, eps)
  DO WHILE  (stackpointer .GE. 1 .OR. active_threads .NE. 0)
! *** ours

! *** Pop next interval off stack.
! *** ours
     active_threads = active_threads + 1
! *** ours
     l   = stack(il,stackpointer)
     r   = stack(ir,stackpointer)
     eps = stack(ie,stackpointer)
     stackpointer = stackpointer - 1

! *** Compute estimates.

     m    = 0.5 * (l + r)
     est1 = 0.5 * (r - l) * (Func(l) + Func(r))
     est2 = 0.5 * ((m - l) * (Func(l) + Func(m)) + &
            (r - m) * (Func(m) + Func(r)))
     abserror = ABS(est2-est1) / 3.0

! *** Check for desired accuracy: push both halves onto the
! *** stack if not accurate enough.

     IF (abserror .LE. eps) THEN
! *** ours
        result = result + est2
        active_threads = active_threads - 1
! *** ours
     ELSE
        IF (stackpointer+2 .GT. STACKSIZE) THEN
           PRINT *, "Stack too small, try STACKSIZE = ", 2*STACKSIZE
           STOP
        END IF

! *** ours
! *** ours
        stackpointer = stackpointer + 1
        stack(il,stackpointer) = l
        stack(ir,stackpointer) = m
        stack(ie,stackpointer) = eps * 0.5

        stackpointer = stackpointer + 1
        stack(il,stackpointer) = m
        stack(ir,stackpointer) = r
        stack(ie,stackpointer) = eps * 0.5
! *** ours
        active_threads = active_threads - 1
! *** ours

     END IF
  END DO
!$omp end parallel
  RETURN
END SUBROUTINE Stackops
       
REAL*8 FUNCTION Func1 (x)

  REAL*8 x, q
  INTEGER n, i

  n = 2000
  q = 1000.0

  DO i=1,n
     q = q - x
  END DO
  if (q.eq.1.0d10) print *,q

  Func1 = x * x
  RETURN
END FUNCTION Func1

REAL*8 FUNCTION Func2(x)

  REAL*8 x, q
  INTEGER n, i

  n = 20000
  q = 1000.0

  DO i=1,n
     q = q - x
  END DO
  if (q.eq.1.0d10) print *,q

  Func2 = x * x
  RETURN
END FUNCTION Func2

REAL*8 FUNCTION Func3 (x)

  REAL*8 x, q
  INTEGER n, i

  n = 200000
  q = 1000.0

  DO i=1,n
     q = q - x
  END DO
  if (q.eq.1.0d10) print *,q

  Func3 = x * x
  RETURN
END FUNCTION Func3
