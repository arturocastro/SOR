PROGRAM Main

! ***  Solution of Laplace's Equation.
! ***
! ***  Uxx + Uyy = 0
! ***  0 <= x <= pi, 0 <= y <= pi
! ***  U(x,pi) = sin(x), U(x,0) = U(0,y) = U(pi,y) = 0
! ***
! ***  then U(x,y) = (sinh(y)*sin(x)) / sinh(pi)
! ***
! ***  Should converge with
! ***            tol = 0.001 and M = 20  in  42 iterations.
! ***   and with tol = 0.001 and M = 100 in 198 iterations.
! *** 

  use timers

  type(TimerType) :: Total

  INTEGER M
  PARAMETER (M = 100)

  REAL*8   PI
  REAL*8   DATAN
  REAL*8   unew(0:M+1,0:M+1), uold(0:M+1,0:M+1)
  REAL*8   solution(0:M+1,0:M+1)
  REAL*8   omega, tol, h
  INTEGER  i, j, iters

  call TimerCreate(Total,"total");

  call TimerOn(Total)

  PI = 4.0D0*DATAN(1.0D0)

  h = PI/FLOAT(M+1)

  DO i=0,M+1
     uold(i,M+1) = SIN(FLOAT(i)*h)
  END DO

  DO i=0,M+1
     DO j=0,M
        uold(i,j) = FLOAT(j)*h*uold(i,M+1)
     END DO
  END DO

  DO i=0,M+1
     DO j=0,M+1
        solution(i,j) = SINH(FLOAT(j)*h)*SIN(FLOAT(i)*h)/SINH(PI)
     END DO
  END DO

  omega = 2.0/(1.0+SIN(PI/FLOAT(M+1)))
  tol   = 0.001

  CALL SOR (unew, uold, solution, omega, tol, m, iters)

  call TimerOff(Total)

  PRINT *, " "
  PRINT *, " Omega = ", omega
  PRINT *, " It took ", iters, " iterations."

  call TimerPrint(Total)

  STOP
END PROGRAM Main


REAL*8 FUNCTION ComputeError (solution, u, m, iters)
  
  INTEGER m, iters
  REAL*8  u(0:m+1,0:m+1)
  REAL*8  solution(0:m+1,0:m+1)

! *** Local variables

  REAL*8  error
  INTEGER i, j

  error = 0.0

  DO j=1,m
     DO i=1,m
        error = MAX(error, ABS(solution(i,j)-u(i,j)))
     END DO
  END DO

!  PRINT *, "On iteration ", iters, " error = ", error

  ComputeError = error

  RETURN
END FUNCTION ComputeError

SUBROUTINE SOR (unew, uold, solution, omega, tol, m, iters)

  INTEGER m, iters
  REAL*8  unew(0:m+1,0:m+1), uold(0:m+1,0:m+1)
  REAL*8  solution(0:m+1,0:m+1)
  REAL*8  omega, tol

! *** Local variables

  REAL*8  error
  INTEGER i, j

! *** External function.

  REAL*8   ComputeError
  EXTERNAL ComputeError

! *** Copy bonudary conditions.

  DO i=0,m+1
     unew(i,m+1) = uold(i,m+1)
     unew(m+1,i) = uold(m+1,i)
     unew(i,  0) = uold(i,  0)
     unew(0,  i) = uold(0,  i)
  END DO

! *** Do SOR until 'tol' satisfied.

  iters = 0
  error = ComputeError (solution, uold, m, iters)

  DO WHILE (error .GE. tol)

! *** Do one iteration of SOR

     !$OMP PARALLEL DO PRIVATE(i, j) SCHEDULE(STATIC,25)
     DO j=1,m
        DO i=1,m
           unew(i,j) = 0.25 * (uold(i - 1, j) + uold(i + 1, j) + uold(i, j - 1) + uold(i, j + 1))
        END DO
     END DO

! *** Copy new to old.

     DO j=1,m
        DO i=1,m
           uold(i,j) = unew(i,j)
        END DO
     END DO

! *** Check error every 20 iterations.

     iters = iters + 1

     IF (MOD(iters,20) .EQ. 0) THEN
        error = ComputeError (solution, uold, m, iters)
     END IF

  END DO

  RETURN
END SUBROUTINE SOR
