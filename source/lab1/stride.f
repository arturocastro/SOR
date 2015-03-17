C
C STRIDE EXAMPLE
C THIS VERSION WRITTEN BY MARK BULL 8/9/97
C
      integer n
      parameter (n=2000,m=20000)
      real a(n,m),junk(2),time1,time2,start
c
C INITIALIZE ARRAY
c
      do j=1,m
         do i=1,n 
            a(i,j) = 0.0 
         end do
      end do
c
      start = etime(junk)
c
C  LOOP NEST 1 
c
      do j=1,m
         do i=1,n 
            a(i,j) = a(i,j) +1.0 
         end do
      end do
c
      time1 = etime(junk) - start
c 
      start = etime(junk)
c
C  LOOP NEST 2
c
      do i=1,n
         do j=1,m 
            a(i,j) = a(i,j) +1.0 
         end do
      end do
c
      time2 = etime(junk) - start
c
      print *, "Time for loop nest 1 = ",time1," seconds"
      print *, "Time for loop nest 2 = ",time2," seconds"
c
      stop
      end
