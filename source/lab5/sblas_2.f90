program sblas
      
  use mpi
  use timers
  type(TimerType) :: vecaddTime,vecsumTime,trimvTime,totalTime

! .. Local Scalars
  integer i,j
! Logical sizes of arrays
  integer n, m
  parameter (n=1024*60,m=512)
  integer k,nsteps
  parameter (nsteps=1000)

! Storage for this processors part of the data.
! (set large enough for single processor to allow runs with 
! different numbers of processors without recompilation. For
! memory efficiency and scalability the space allocated should
! be around 1/pth of the total memory requirement)

  integer ns, ms
  parameter (ns=1024*60,ms=512)

  real*8 a(ns),b(ns),x(ns),y(ms),l(ms,ms), a_g(ns)
  common /heap/ a,b,x,y,l,a_g
  real*8 res, glob_res

! MPI specific data

  integer myid, numprocs, rc, ierr
! Local array sizes (on a particular processor)
  integer nx, nnx, mx, nmx, modx, pdn, pdm

! ..Intrinsic functions
  intrinsic MOD

! .. Executable statements

! MPI initialisation

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  print *, "Process ", myid, " of ", numprocs, " is alive"

  call TimerCreate(vecaddTime,"vecadd")
  call TimerCreate(vecsumTime,"vecsum")
  call TimerCreate(trimvTime,"trimv")
  if (myid.eq.0) call TimerCreate(totalTime,"total")

! Loop over whole program to get reasonable execution time
! (simulates a time-steping code)

! start total timer on master (proc 0)

  if (myid.eq.0) call TimerOn(totalTime)

  do k=1,nsteps
! start total timer on master (proc 0)

! Work out size of data on this processor and number of processors
! with data

     nx = (n+numprocs-1)/numprocs
     nnx = nx
     pdn = (n+nx-1)/nx
     if (myid.eq.pdn-1) then
        modx = MOD(n,nx)
        if (modx.ne.0) nx = modx
     else if (myid.ge.pdn) then
        nx = 0
     endif

     mx = (m+numprocs-1)/numprocs
     nmx = mx
     pdm = (m+mx-1)/mx
     if (myid.eq.pdm-1) then
        modx = MOD(m,mx)
        if (modx.ne.0) mx = modx
     else if (myid.ge.pdm) then
        mx = 0
     endif

! 'Update' my data

     call vector_update(nx,x,a,b)

! vector addition operation

     call TimerOn(vecaddTime)
     call vecadd(nx,x,a,b)
     call TimerOff(vecaddTime)

! sum operation

     call TimerOn(vecsumTime)
     call vecsum(nx,res,x)
     call TimerOff(vecsumTime)

! Now calculate the sum of the partial results 

     glob_res = 0
     call MPI_ALLREDUCE(res, glob_res, 1, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ierr)

! 'Update' lower triangular matrix and y vector

     call matrix_update(ms,y,l,m)

! triangular matrix times vector operation

! perform the mv operation on the local data.
! Note, array a was distributed block, but now all processors require
! access to it. An all-to-all gather data exchange copies the data 
! to everyone.

! Get a copy of a for each processor.

     call MPI_ALLGather(a, nnx, MPI_DOUBLE_PRECISION, a_g, nnx, &
          MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, &
          ierr)

! Perform the matrix vector operation on the local data
! Pass in the dimension of the local data and the global
! size of the array (for interleaved scheduling)

     call TimerOn(trimvTime)
     call trimv (ms,y,l,a_g,m)
     call TimerOff(trimvTime)

! end time-step loop
  end do

! print out some results

  if (myid .eq. 0) then
     call TimerOff(totalTime)
     print *,'x(1) = ',x(1)
     print *,'Sum of x() = ',glob_res
     print *,'y(1) = ',y(1)
     call TimerPrint(totalTime)
     
  end if

  print *,'x(n) = ', x(nx), 'myid = ', myid
  print *,'y(m) = ', y(mx), 'myid = ', myid

  print *,'myid= ',myid
  call TimerPrint(vecaddTime)
  call TimerPrint(vecsumTime)
  call TimerPrint(trimvTime)

  call MPI_FINALIZE(ierr)
end program sblas

! Subroutines

subroutine vector_update (n,x,a,b)

! .. Scalar arguments
  integer n

! .. Array arguments
  real*8 x(n),a(n),b(n)

! .. Local scalars
  integer i

! .. External functions

  do i=1,n
     x(i) = 0.
     a(i) = 3.142
     b(i) = 3.142
  end do

  return
end subroutine vector_update

subroutine matrix_update (ms,y,l,m)

  USE MPI

! .. Scalar arguments
  integer ms,m

! .. Array arguments
  real*8 y(ms),l(ms,ms)

! .. Local scalars
  integer i,j,ii,ierr
  integer myid,numprocs

! .. External functions

  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! Use interleaved distribution (for load balance). Each processor
! has m/p rows (max) of l and the corresponding elements of y

  ii=1
  do i=myid+1,m,numprocs
     y(ii) = 0.
     do j=1,i
        l(j,ii) = 2.
     end do
     ii=ii+1
  end do

  return
end subroutine matrix_update

subroutine vecadd (n,x,a,b)

! .. Scalar arguments
  integer n

! .. Array arguments
  real*8 x(n),a(n),b(n)

! .. Local scalars
  integer i

! .. External functions
  
  do i=1,n
     x(i) = a(i) + b(i)
  end do
  return
end subroutine vecadd

subroutine vecsum (n,sum,a)

! .. Scalar arguments 
  integer n
  real*8 sum

! .. Array arguments
  real*8 a(n)

! .. Local scalars
  integer i

  sum = 0.0
  do i=1,n
     sum = sum +a(i)
  end do
  return
end subroutine vecsum

subroutine trimv (ms,y,l,a,m)

  USE MPI

! .. Scalar arguments
  integer ms,m

! .. Array arguments
  real*8 y(ms),l(ms,ms),a(ms)

! .. Local scalars
  integer i,j,ii,ierr
  integer myid,numprocs

! .. External functions
  
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

  ii=0
  do i=myid+1,m,numprocs
     ii=ii+1
     do j=1,i
        y(ii) = y(ii) + l(j,ii)*a(j)
     end do
  end do

  return
end subroutine trimv
