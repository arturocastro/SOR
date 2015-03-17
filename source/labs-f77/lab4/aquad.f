      PROGRAM TestQuad

      use timers
      type(TimerType) :: func1Time,func2Time,func3Time

      REAL*8 left, right, result, eps, Quad
      REAL*8     Func1, Func2, Func3
      EXTERNAL   Func1, Func2, Func3
C
C     Initialise timers
      call TimerCreate(func1Time,"Func1")
      call TimerCreate(func2Time,"Func2")
      call TimerCreate(func3Time,"Func3")
C
      call TimerOn(func1Time)
      left  = 0.0
      right = 1.0
      eps   = 1.0E-7
      result = Quad (Func1, left, right, eps)
      call TimerOff(func1Time)
C
      PRINT *, "Result for Func1 = ", result
      call TimerPrint(func1Time)
C
      call TimerOn(func2Time)
      left  = 0.0
      right = 1.0
      eps   = 1.0E-7
      result = Quad (Func2, left, right, eps)
      call TimerOff(func2Time)
C
      PRINT *, "Result for Func2 = ", result
      call TimerPrint(func2Time)
C
      call TimerOn(func3Time)
      left  = 0.0
      right = 1.0
      eps   = 1.0E-7
      result = Quad (Func3, left, right, eps)
      call TimerOff(func3Time)
C
      PRINT *, "Result for Func3 = ", result
      call TimerPrint(func3Time)
C
      END
C
      REAL*8 FUNCTION Quad (Func, left, right, epsilon)
      REAL*8 left, right, epsilon, Func
C
C *** Adaptive Quadrture
C
      INTEGER STACKSIZE,IL,IR,IE
      PARAMETER (STACKSIZE = 1000,IL=1,IR=2,IE=3)

      REAL*8 stack  (3,STACKSIZE)
      INTEGER stackpointer
C
      REAL*8 result, abserror, l, r, m, est1, est2, eps
C
      stackpointer = 1
      stack(il,stackpointer) = left
      stack(ir,stackpointer) = right
      stack(ie,stackpointer) = epsilon
C
      result = 0.0
C
      call Stackops(stackpointer,stack,result,Func) 
C
      Quad = result
C
      RETURN
      END

      SUBROUTINE Stackops(stackpointer,stack,result,Func)
C
      INTEGER STACKSIZE,IL,IR,IE
      PARAMETER (STACKSIZE = 1000,IL=1,IR=2,IE=3)

      REAL*8 stack(3,STACKSIZE)
      INTEGER stackpointer
C
      REAL*8 result, Func
      REAL*8 abserror, l, r, m, est1, est2, eps
C
      DO WHILE (stackpointer .GE. 1)
C
C *** Pop next interval off stack.
C
         l   = stack(il,stackpointer)
         r   = stack(ir,stackpointer)
         eps = stack(ie,stackpointer)
         stackpointer = stackpointer - 1
C
C *** Compute estimates.
C
         m    = 0.5 * (l + r)
         est1 = 0.5 * (r - l) * (Func(l) + Func(r))
         est2 = 0.5 * ((m - l) * (Func(l) + Func(m)) + 
     1                 (r - m) * (Func(m) + Func(r)))
         abserror = ABS(est2-est1) / 3.0
C
C *** Check for desired accuracy: push both halves onto the
C *** stack if not accurate enough.
C
         IF (abserror .LE. eps) THEN
            result = result + est2
         ELSE
            IF (stackpointer+2 .GT. STACKSIZE) THEN
               PRINT *, "Stack too small, try STACKSIZE = ", 2*STACKSIZE
               STOP
            END IF
C
            stackpointer = stackpointer + 1
            stack(il,stackpointer) = l
            stack(ir,stackpointer) = m
            stack(ie,stackpointer) = eps * 0.5
C
            stackpointer = stackpointer + 1
            stack(il,stackpointer) = m
            stack(ir,stackpointer) = r
            stack(ie,stackpointer) = eps * 0.5
         END IF
      END DO
C
      RETURN
      END
C       
      REAL*8 FUNCTION Func1 (x)
C
      REAL*8 x, q
      INTEGER n, i
C
      n = 2000
      q = 1000.0
C
      DO i=1,n
         q = q - x
      END DO
      if (q.eq.1.0d10) print *,q
C
      Func1 = x * x
      RETURN
      END
C
      REAL*8 FUNCTION Func2(x)
C
      REAL*8 x, q
      INTEGER n, i
C
      n = 20000
      q = 1000.0
C
      DO i=1,n
         q = q - x
      END DO
      if (q.eq.1.0d10) print *,q
C
      Func2 = x * x
      RETURN
      END
C
      REAL*8 FUNCTION Func3 (x)
C
      REAL*8 x, q
      INTEGER n, i
C
      n = 200000
      q = 1000.0
C
      DO i=1,n
         q = q - x
      END DO
      if (q.eq.1.0d10) print *,q
C
      Func3 = x * x
      RETURN
      END
