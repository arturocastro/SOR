

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <string.h>
    #include <sys/time.h>
    #include <mpi.h>
     
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
     
#define TOLERANCE 0.00002        /* termination criterion */
#define MAGIC 0.8                /* magic factor */

#define M 20
     
int N; /* problem size */
double **G; /* the grid */
double **buff;
MPI_Status istatus;
MPI_Request ireqs1, ireqs2, ireqr1, ireqr2;
int mynode, totalnodes;
int jsta, jend, inext, iprev;
double *bufs1, *bufr1;
double *bufs2, *bufr2;
double stopdiff, maxdiff, diff; /* differences btw grid points in iters */



void shift(const int iflg)
{
    double bufs1a[M + 2] = {0};
    double bufs2a[M + 2] = {0};
    
    double bufr1a[M + 2] = {0};
    double bufr2a[M + 2] = {0};

    const int is1 = ((jsta + iflg) % 2) + 1;
    const int is2 = ((jend + iflg) % 2) + 1;

    const int ir1 = fabs(3 - is1);
    const int ir2 = fabs(3 - is2);

    int i, icnt1=0, icnt2=0;

    if(mynode != 0)
    {
	icnt1 = 0;

	for(i = is1; i < M + 2; i += 2)
	{
	    bufs1a[icnt1] = G[i][jsta];
	    ++icnt1;
	}
    }

    if (mynode != totalnodes - 1)
    {
	icnt2 = 0;

#ifdef DEBUG
	printf("rank = %d, ied=%d\n", rank, jend);
#endif
	
	for(i = is2; i < M + 2; i += 2)
	{
	    bufs2a[icnt2] = G[i][jend];

#ifdef DEBUG
	    printf("s%f ", bufs2[icnt2]);
#endif
	    ++icnt2;
	}
    }

    //if(iflg == 0)
	//printf("iflg=%d\nis1 = %d, is2 = %d\nir1 = %d, ir2=%d\n", iflg, is1, is2, ir1, ir2);

    MPI_Isend(bufs1a, icnt1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqs1);
    MPI_Isend(bufs2a, icnt2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqs2);
    
    MPI_Irecv(bufr1a, M + 2, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqr1);
    MPI_Irecv(bufr2a, M + 2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqr2);
    
    MPI_Wait(&ireqs1, &istatus);
    MPI_Wait(&ireqs2, &istatus);
    
    MPI_Wait(&ireqr1, &istatus);
    MPI_Wait(&ireqr2, &istatus);

    int icnt;
   
    if(mynode != 0)
    {
	icnt = 0;

#ifdef DEBUF
	printf("rank = %d, jsta=%d\n", rank, jsta);
#endif

	for(i = ir1; i < M + 2; i += 2)
	{
	    G[i][jsta - 1] = bufr1a[icnt];
	    
#ifdef DEBUG
	    printf("r%f ", bufr1a[icnt]);
#endif
	    ++icnt;
	}
    }

    if(mynode != totalnodes - 1)
    {
	icnt = 0;

	for(i = ir2; i < M + 2; i += 2)
	{
	    G[i][jend + 1] = bufr2a[icnt];
	    ++icnt;
	}
    }
}
     
int even (int i)
{
    return !(i & 1);
}
     
double stencil (double** G, int row, int col)
{
    return (G[row-1][col] + G[row+1][col] + G[row][col-1] + G[row][col+1] ) / 4.0;
}
     
void alloc_grid(double ***Gptr, int N)
{
    int i;
    double** G = (double**)malloc(N*sizeof(double*));
    if ( G == 0 ) {
	fprintf(stderr, "malloc failed\n");
	exit(42);
    }
     
    for (i = 0; i<N; i++) { /* malloc the own range plus one more line */
	/* of overlap on each side */
	G[i] = (double*)malloc(N*sizeof(double));
	if ( G[i] == 0 ) {
	    fprintf(stderr, "malloc failed\n");
	    exit(42);
	}
    }
     
    *Gptr = G;
}
     
void init_grid(double **G, int N)
{
    int i, j;
     
    /* initialize the grid */
    for (i = 0; i < N; i++){
	for (j = 0; j < N; j++){
	    if (i == 0) G[i][j] = 4.56;
	    else if (i == N-1) G[i][j] = 9.85;
	    else if (j == 0) G[i][j] = 7.32;
	    else if (j == N-1) G[i][j] = 6.88;
	    else G[i][j] = 0.0;
	}
    }
}
     
void print_grid(double** G, int N)
{
    int i, j;
     
    for ( i = 1 ; i < N-1 ; i++ ) {
	for ( j = 1 ; j < N-1 ; j++ ) {
	    printf("%10.3f ", G[i][j]);
	}
	printf("\n");
    }
}
     
void range(int n1, int n2, int nprocs, int irank, int *ista, int *iend) {
    float iwork1 = (n2 - n1 + 1) / nprocs;
    float iwork2 = (n2 - n1 + 1) % nprocs;
    int start, end;
     
    start = (int) (irank * iwork1 + n1 + MIN(irank, iwork2));
    end = (int) (start + iwork1 - 1);
     
    if (iwork2 > irank) {
	end = end + 1;
    }
     
    *ista = start;
    *iend = end;
}
     
void exchange(int phase) {
    int is1 = ((jsta + phase) % 2) + 1;
    int is2 = ((jend + phase) % 2) + 1;
    int ir1 = fabs(3 - is1);
    int ir2 = fabs(3 - is2);
    int icnt = 0;
    int icnt1 = 0, icnt2 = 0;
    int m = N - 1;
    int i;
     
    if (mynode != 0) {
	icnt1 = 0;
	for (i = is1; i <= m; i += 2) {
	    bufs1[icnt1] = G[i][jsta];
	    icnt1 = icnt1 + 1;
	}
    }
     
    if (mynode != (totalnodes - 1)) {
	icnt2 = 0;
	for (i = is2; i <= m; i += 2) {
	    bufs2[icnt2] = G[i][jend];
	    icnt2 = icnt2 + 1;
	}
    }
     
    MPI_Isend(bufs1, icnt1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqs1);
    MPI_Isend(bufs2, icnt2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqs2);
     
    MPI_Irecv(bufr1, N, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqr1);
    MPI_Irecv(bufr2, N, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqr2);
     
    MPI_Wait(&ireqs1, &istatus);
    MPI_Wait(&ireqs2, &istatus);
    MPI_Wait(&ireqr1, &istatus);
    MPI_Wait(&ireqr2, &istatus);
     
    if (mynode != 0) {
	icnt = 0;
	for (i = ir1; i <= m; i += 2) {
	    G[i][jsta - 1] = bufr1[icnt];
	    icnt = icnt + 1;
	}
    }
     
    if (mynode != (totalnodes - 1)) {
	icnt = 0;
	for (i = ir2; i <= m; i += 2) {
	    G[i][jend + 1] = bufr2[icnt];
	    icnt = icnt + 1;
	}
    }
     
}
     
int main (int argc, char *argv[])
{
    int ncol,nrow;             /* number of rows and columns */
    double Gnew,r,omega,       /* some float values */
	stopdiff,maxdiff,diff;  /* differences btw grid points in iters */
    int i,j,phase,iteration;   /* counters */
    int print = 0;
    double inittime, totaltime;
     
    /* Initializing MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
     
    /* set up problem size */
    N = 20;
     
    for(i=1; i<argc; i++) {
	if(!strcmp(argv[i], "-print")) {
	    print = 1;
	} else {
	    N = atoi(argv[i]);
	}
    }
     
    if(mynode == 0) {
	fprintf(stderr, "Running %d x %d SOR\n", N, N);
    }
     
    N += 2; /* add the two border lines */
            /* finally N*N (from argv) array points will be computed */
     
            /* set up a quadratic grid */
    ncol = nrow = N;
    r        = 0.5 * ( cos( M_PI / ncol ) + cos( M_PI / nrow ) );
    omega    = 2.0 / ( 1 + sqrt( 1 - r * r ) );
    stopdiff = TOLERANCE / ( 2.0 - omega );
    omega   *= MAGIC;
     
    alloc_grid(&G, N);
    alloc_grid(&buff, N);
    init_grid(G, N);
     
    // send receive buffers for left neighbor processes (iprev)
    bufs1 = (double *) malloc(N * sizeof (double));
    bufr1 = (double *) malloc(N * sizeof (double));
     
    // send receive buffers for right neighbor processes (inext)
    bufs2 = (double *) malloc(N * sizeof (double));
    bufr2 = (double *) malloc(N * sizeof (double));
     
    // start - end matrix array cols for each processes
    range(1, N - 1, totalnodes, mynode, &jsta, &jend);
    
    printf("numtasks = %d\nrank %d: from %d to %d\n", totalnodes, mynode, jsta, jend);
     
    inext = mynode + 1;
    iprev = mynode - 1;
     
    if (inext == totalnodes)
	inext = MPI_PROC_NULL;
     
    if (iprev == -1)
	iprev = MPI_PROC_NULL;
     
    inittime = MPI_Wtime();
     
    /* now do the "real" computation */
    iteration = 0;
    do {
	maxdiff = 0.0;
	for ( phase = 0; phase < 2 ; phase++){
     
	    shift(phase);
     
	    for ( i = 1 ; i < N-1 ; i++ ){
		for ( j = jsta + (even(i) ^ phase); j <= jend ; j += 2 ){
		    Gnew = stencil(G,i,j);
		    diff = fabs(Gnew - G[i][j]);
		    if ( diff > maxdiff )
			maxdiff = diff;
		    G[i][j] = G[i][j] + omega * (Gnew-G[i][j]);
		}
	    }
	}
     
	MPI_Allreduce(&maxdiff, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	maxdiff = diff;
     
	iteration++;
    } while (maxdiff > stopdiff);
     
    totaltime = MPI_Wtime() - inittime;
     
    MPI_Barrier(MPI_COMM_WORLD);
     
    if(print == 1) {
	printf("Node: %d\n", mynode);
	print_grid(G, N);
	printf("\n");
    }
     
    MPI_Finalize();
     
printf("%f", totaltime);
    return 0;
}

