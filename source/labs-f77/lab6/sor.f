      PROGRAM Main
C
C ***  Solution of Laplace's Equation.
C ***
C ***  Uxx + Uyy = 0
C ***  0 <= x <= pi, 0 <= y <= pi
C ***  U(x,pi) = sin(x), U(x,0) = U(0,y) = U(pi,y) = 0
C ***
C ***  then U(x,y) = (sinh(y)*sin(x)) / sinh(pi)
C ***
C ***  Should converge with
C ***            tol = 0.001 and M = 20  in  42 iterations.
C ***   and with tol = 0.001 and M = 100 in 198 iterations.
C *** 
C
      INTEGER M
      PARAMETER (M   = 100)
C
      REAL*8   PI
      REAL*8   DATAN
      REAL*8   unew(0:M+1,0:M+1), uold(0:M+1,0:M+1)
      REAL*8   solution(0:M+1,0:M+1)
      REAL*8   omega, tol, h
      INTEGER  i, j, iters
C
      PI = 4.0D0*DATAN(1.0D0)
C
      h = PI/FLOAT(M+1)
C
      DO i=0,M+1
         uold(i,M+1) = SIN(FLOAT(i)*h)
      END DO
C
      DO i=0,M+1
         DO j=0,M
            uold(i,j) = FLOAT(j)*h*uold(i,M+1)
         END DO
      END DO
C
      DO i=0,M+1
         DO j=0,M+1
            solution(i,j) = SINH(FLOAT(j)*h)*SIN(FLOAT(i)*h)/SINH(PI)
         END DO
      END DO
C
      omega = 2.0/(1.0+SIN(PI/FLOAT(M+1)))
      tol   = 0.001
C
      CALL SOR (unew, uold, solution, omega, tol, m, iters)
C
      PRINT *, " "
      PRINT *, " Omega = ", omega
      PRINT *, " It took ", iters, " iterations."
C
      STOP
      END
C

      REAL*8 FUNCTION ComputeError (solution, u, m, iters)
C
      INTEGER m, iters
      REAL*8  u(0:m+1,0:m+1)
      REAL*8  solution(0:m+1,0:m+1)
C
C *** Local variables
C
      REAL*8  error
      INTEGER i, j

      error = 0.0
C
      DO j=1,m
         DO i=1,m
            error = MAX(error, ABS(solution(i,j)-u(i,j)))
         END DO
      END DO
C
      PRINT *, "On iteration ", iters, " error = ", error
C
      ComputeError = error
C
      RETURN
      END
C
      SUBROUTINE SOR (unew, uold, solution, omega, tol, m, iters)
C
      INTEGER m, iters
      REAL*8  unew(0:m+1,0:m+1), uold(0:m+1,0:m+1)
      REAL*8  solution(0:m+1,0:m+1)
      REAL*8  omega, tol
C
C *** Local variables
C
      REAL*8  error
      INTEGER i, j
C
C *** External function.
C
      REAL*8   ComputeError
      EXTERNAL ComputeError
C
C *** Copy bonudary conditions.
C
      DO i=0,m+1
         unew(i,m+1) = uold(i,m+1)
         unew(m+1,i) = uold(m+1,i)
         unew(i,  0) = uold(i,  0)
         unew(0,  i) = uold(0,  i)
      END DO
C
C *** Do SOR until 'tol' satisfied.
C
      iters = 0
      error = ComputeError (solution, uold, m, iters)
C
      DO WHILE (error .GE. tol)
C
C *** Do one iteration of SOR
C
         DO j=1,m
            DO i=1,m
               unew(i,j) = uold(i,j) + 0.25*omega*
     1                    (unew(i-1,j) + unew(i,j-1) +
     2                     uold(i+1,j) + uold(i,j+1) -
     3                     4.0*uold(i,j))
            END DO
         END DO
C
C *** Copy new to old.
C
         DO j=1,m
            DO i=1,m
               uold(i,j) = unew(i,j)
            END DO
         END DO
C
C *** Check error every 20 iterations.
C
         iters = iters + 1
C
         IF (MOD(iters,20) .EQ. 0) THEN
            error = ComputeError (solution, uold, m, iters)
         END IF
C
      END DO
C
      RETURN
      END
