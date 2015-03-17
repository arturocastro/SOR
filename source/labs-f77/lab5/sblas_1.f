      program sblas
c    
      use mpi
      use timers
      type(TimerType) :: vecaddTime,vecsumTime,trimvTime,totalTime

      integer i,j,n,m,k,nsteps
      integer nx
      integer status(MPI_STATUS_SIZE)
      parameter (n=1024*60,m=512)
      parameter (nsteps=1000)
      real*8 a(n),b(n),x(n),y(m),l(m,m)
      real*8 buff(n)
      common /heap/ a,b,x,y,l
      real*8 res
c
c MPI specific data
c
      integer myid, numprocs, rc, ierr
c
c ..Intrinsic functions
c
c .. Executable statements
c
c MPI initialisation
c
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      print *, "Process ", myid, " of ", numprocs, " is alive"
c
      call TimerCreate(vecaddTime,"vecadd")
      call TimerCreate(vecsumTime,"vecsum")
      call TimerCreate(trimvTime,"trimv")
      if (myid.eq.0) call TimerCreate(totalTime,"total")
c
c Loop over whole program to get reasonable execution time
c (simulates a time-steping code)
c
c start total timer on master (proc 0)
c
      if (myid.eq.0) call TimerOn(totalTime)
c
      do k=1,nsteps
c
c 'Update' vector data on processor 0
c
      if (myid.eq.0) call vector_update(n,x,a,b)
c
c Work out size of data on this processor (not a general solution!)
c
      nx = n/numprocs
c
c Distribute the vectors a and b to other processors
c
      if (myid.eq.0) then
c MASTER
c send parts of the vectors a and b to other processors (x not reqd)
c
         do i=1,numprocs-1
            call MPI_SEND(a(i*nx+1), nx, MPI_DOUBLE_PRECISION, 
     +           I, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(b(i*nx+1), nx, MPI_DOUBLE_PRECISION, 
     +           I, 0, MPI_COMM_WORLD, ierr)
         end do
      else
c SLAVE
c receive data
c
         call MPI_RECV(a, nx, MPI_DOUBLE_PRECISION,
     +        0, 0, MPI_COMM_WORLD,
     +        status, ierr)
         call MPI_RECV(b, nx, MPI_DOUBLE_PRECISION,
     +        0, 0, MPI_COMM_WORLD,
     +        status, ierr)
      end if 
c
c vector addition operation: each processor adds its own parts 
c Note: only master time is collected...
c
      call TimerOn(vecaddTime)
      call vecadd (nx,x,a,b)
      call TimerOff(vecaddTime)
c
c Gather result vector on processor 0
c
      if (myid.eq.0) then
c MASTER
c receive parts of result into the appropriate parts of my x
c
         do i=1,numprocs-1
            call MPI_RECV(x(i*nx+1), nx, MPI_DOUBLE_PRECISION,
     +           i, 0, MPI_COMM_WORLD,
     +           status, ierr)
         end do
c
      else
c SLAVE
c send my part of result vector
c
         call MPI_SEND(x, nx, MPI_DOUBLE_PRECISION, 
     +        0, 0, MPI_COMM_WORLD, ierr)
c
      end if 
c
c Rest of the code is performed on processor 0
c
      if (myid.eq.0) then
c
c sum operation
c
      call TimerOn(vecsumTime)
      call vecsum(n,res,x)
      call TimerOff(vecsumTime)
c
c 'Update' lower triangular matrix and y vector
c
      call matrix_update(m,y,l)
c
c triangular matrix times vector operation
c
      call TimerOn(trimvTime)
      call trimv (m,y,l,a)
      call TimerOff(trimvTime)
c
      end if
c
c end time-step loop
      end do
c
c print out some results
c
      if (myid.eq.0) then
         call TimerOff(totalTime)
         print *,'x(1) = ',x(1),'x(n) = ', x(n)
         print *,'Sum of x() = ',res
         print *,'y(1) = ',y(1),'y(m) = ', y(m)
         call TimerPrint(vecsumTime)
         call TimerPrint(trimvTime)
         call TimerPrint(totalTime)
c
c end of processor 0 processing
      end if
c
      print *,'myid= ', myid
      call TimerPrint(vecaddTime)
c
      call MPI_FINALIZE(ierr)
      end
c
c Subroutines
c
      subroutine vector_update(n,x,a,b)
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
      subroutine matrix_update(m,y,l)
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
      subroutine vecsum(n,sum,a)
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
      do 10 i=1,n
        sum = sum +a(i)
10    continue
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

