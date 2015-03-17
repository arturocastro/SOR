      program sblas
c      
      use mpi
      use timers
      type(TimerType) :: vecaddTime,vecsumTime,trimvTime,totalTime
c
c .. Local Scalars
      integer i,j
c Logical sizes of arrays
      integer n, m
      parameter (n=1024*60,m=512)
      integer k,nsteps
      parameter (nsteps=1000)
c
c Storage for this processors part of the data.
c (set large enough for single processor to allow runs with 
c different numbers of processors without recompilation. For
c memory efficiency and scalability the space allocated should
c be around 1/pth of the total memory requirement)
c
      integer ns, ms
      parameter (ns=1024*60,ms=512)
c 
      real*8 a(ns),b(ns),x(ns),y(ms),l(ms,ms), a_g(ns)
      common /heap/ a,b,x,y,l,a_g
      real*8 res, glob_res
c
c MPI specific data
c
      integer myid, numprocs, rc, ierr
c Local array sizes (on a particular processor)
      integer nx, nnx, mx, nmx, modx, pdn, pdm
c
c ..Intrinsic functions
      intrinsic MOD
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
c
c Loop over whole program to get reasonable execution time
c (simulates a time-steping code)
c
c start total timer on master (proc 0)
c
      if (myid.eq.0) call TimerOn(totalTime)
c
      do k=1,nsteps
c start total timer on master (proc 0)
c
c Work out size of data on this processor and number of processors
c with data
c
      nx = (n+numprocs-1)/numprocs
      nnx = nx
      pdn = (n+nx-1)/nx
      if (myid.eq.pdn-1) then
         modx = MOD(n,nx)
         if (modx.ne.0) nx = modx
      else if (myid.ge.pdn) then
         nx = 0
      endif
c
      mx = (m+numprocs-1)/numprocs
      nmx = mx
      pdm = (m+mx-1)/mx
      if (myid.eq.pdm-1) then
         modx = MOD(m,mx)
         if (modx.ne.0) mx = modx
      else if (myid.ge.pdm) then
         mx = 0
      endif
c
c 'Update' my data
c
      call vector_update(nx,x,a,b)
c
c vector addition operation
c
      call TimerOn(vecaddTime)
      call vecadd(nx,x,a,b)
      call TimerOff(vecaddTime)
c
c sum operation
c
      call TimerOn(vecsumTime)
      call vecsum(nx,res,x)
      call TimerOff(vecsumTime)
c
c Now calculate the sum of the partial results 
c
      glob_res = 0
      call MPI_ALLREDUCE(res, glob_res, 1, MPI_DOUBLE_PRECISION,
     +                   MPI_SUM, MPI_COMM_WORLD, ierr)
c
c 'Update' lower triangular matrix and y vector
c
      call matrix_update(ms,y,l,m)
c
c triangular matrix times vector operation
c
c perform the mv operation on the local data.
c Note, array a was distributed block, but now all processors require
c access to it. An all-to-all gather data exchange copies the data 
c to everyone.
c
c Get a copy of a for each processor.
c
      call MPI_ALLGather(a, nnx, MPI_DOUBLE_PRECISION, a_g, nnx,
     +                   MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,
     +                   ierr)
c
c Perform the matrix vector operation on the local data
c Pass in the dimension of the local data and the global
c size of the array (for interleaved scheduling)
c
      call TimerOn(trimvTime)
      call trimv (ms,y,l,a_g,m)
      call TimerOff(trimvTime)
c
c end time-step loop
	end do
c
c print out some results
c
      if (myid .eq. 0) then
         call TimerOff(totalTime)
         print *,'x(1) = ',x(1)
         print *,'Sum of x() = ',glob_res
         print *,'y(1) = ',y(1)
         call TimerPrint(totalTime)

      end if
c
      print *,'x(n) = ', x(nx), 'myid = ', myid
      print *,'y(m) = ', y(mx), 'myid = ', myid
c
      print *,'myid= ',myid
      call TimerPrint(vecaddTime)
      call TimerPrint(vecsumTime)
      call TimerPrint(trimvTime)

      call MPI_FINALIZE(ierr)
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
      subroutine matrix_update (ms,y,l,m)
c
      include 'mpif.h'
c
c .. Scalar arguments
      integer ms,m
c
c .. Array arguments
      real*8 y(ms),l(ms,ms)
c
c .. Local scalars
      integer i,j,ii,ierr
      integer myid,numprocs
c
c .. External functions
c
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c
c Use interleaved distribution (for load balance). Each processor
c has m/p rows (max) of l and the corresponding elements of y
c
      ii=1
      do i=myid+1,m,numprocs
         y(ii) = 0.
         do j=1,i
            l(j,ii) = 2.
         end do
         ii=ii+1
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
      subroutine trimv (ms,y,l,a,m)
c
      include 'mpif.h'
c
c .. Scalar arguments
      integer ms,m
c
c .. Array arguments
      real*8 y(ms),l(ms,ms),a(ms)
c
c .. Local scalars
      integer i,j,ii,ierr
      integer myid,numprocs
c
c .. External functions
c
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c
      ii=0
      do i=myid+1,m,numprocs
         ii=ii+1
         do j=1,i
            y(ii) = y(ii) + l(j,ii)*a(j)
         end do
      end do
c
      return
      end
