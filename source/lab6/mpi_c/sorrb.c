#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
 
#define M 500

int numtasks, rank, inext, iprev, jsta, jend;

//double solution[M + 2][M + 2] = {{ 0 }};

double uold[M + 2][M + 2] = {{ 0 }};

MPI_Request ireqs1;
MPI_Request ireqs2;

MPI_Request ireqr1;
MPI_Request ireqr2;

MPI_Status istatus;

int sor(double uold[][M + 2], const double omega, const double stopdiff, const int m);
void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict jsta, int * restrict jend);
void shift(const int iflg);
inline int min(const int a, const int b);

inline int even(const int i)
{
    return !(i & 1);
}

void printmat()
{
    int i, j;
    for(i = 0; i < M + 2; ++i)
    {
	for(j = 0; j < M + 2; ++j)
	{
	    printf("%.2f ", uold[i][j]);
	}
	puts("\n");
    }
}

void printmatrank()
{
    int i, j;
    for (i = jsta; i <= jend; ++i)
    {
	for (j = 0; j < M + 2; ++j)
	{
	    printf("%.2f ", uold[i][j]);
	}
	puts("\n");
    }
}

int main(int argc, char * argv[])
{
    const int m = M;

    int ierr = MPI_Init(&argc, &argv);

    if(ierr != MPI_SUCCESS)
    {
	perror("MPI init failed. Terminating T800.\n");
	exit(1);
    }

    int i, j;

    const double begin = MPI_Wtime();

    const double pi = 4.0 * atan(1.0);

    const double h = pi / (m + 1);

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(i = 0; i < m + 2; ++i)
    {
       uold[i][M + 1] = sin(i * h);
    }

    for(i = 0; i < m + 2; ++i)
    {
	for(j = 0; j < m + 1; ++j)
        {
            uold[i][j] = j * h * uold[i][M + 1];
	    //uold[i][j] = (i + 1) * 100 + j;
        }
    }

    /*for(i = 0; i < m + 2; ++i)
    {
        for(j = 0; j < m + 2; ++j)
        {
            solution[i][j] = sinh(j * h) * sin(i * h) / sinh(pi);
        }
    }*/

    para_range(1, m + 1, numtasks, rank, &jsta, &jend);

    printf("numtasks = %d\nrank %d: from %d to %d\n", numtasks, rank, jsta, jend);

    inext = rank + 1;
    iprev = rank - 1;

    if(inext == numtasks)
    {
	inext = MPI_PROC_NULL;
    }

    if(iprev == -1)
    {
	iprev = MPI_PROC_NULL;
    }

    const double omega = 2.0 / ( 1.0 + sin(pi / (m + 1)) );
    const double tol = 0.001;

    const double stopdiff = tol / (2.0 - omega);

    const int iters = sor(uold, omega, stopdiff, m);

    const double end = MPI_Wtime();

    MPI_Finalize();

    if(rank == 0)
    {
	printf(" \n");
	printf(" Omega = %f\n", omega);
	printf(" It took %d iterations.\n", iters);
	
	printf("Total time = %f\n\n\n", end - begin);
    }

    return 0;
}

int sor(double uold[][M + 2], const double omega, const double stopdiff, const int m)
{
    int i, j;

    int iters = 0;
    double temp;
    double maxdiff = 0.0;
    double diff;

    // Do SOR until 'tol' is satisfied
    do
    {
	maxdiff = 0.0;
#ifdef DEBUG
	if(rank == 1)
	{
	    printf("jsta = %d, jend = %d\n", jsta, jend);

	    printmat();	    
	}
#endif

	// Both for loops account to one iteration of SOR

	/* Update only Red
	  +----+----+----+----+----+----+
	  | B  | R* | B  | R* | B  | R* |
	  +----+----+----+----+----+----+
	  | R* | B  | R* | B  | R* | B  |
	  +----+----+----+----+----+----+
	  | B  | R* | B  | R* | B  | R* |
	  +----+----+----+----+----+----+
	  | R* | B  | R* | B  | R* | B  |
	  +----+----+----+----+----+----+
	  | B  | R* | B  | R* | B  | R* |
	  +----+----+----+----+----+----+
	*/
	shift(0);

        for(i = 1; i < m + 1; ++i)
        {
            for(j = jsta + (even(i) ^ 0); j <= jend; j += 2)
            {
 		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		diff = fabs(temp - uold[i][j]);
		if (diff > maxdiff) { maxdiff = diff; }
		uold[i][j] += omega * temp;
            }
        }

	/* Update only Black
	  +----+----+----+----+----+----+
	  | B* | R* | B* | R* | B* | R* |
	  +----+----+----+----+----+----+
	  | R* | B* | R* | B* | R* | B* |
	  +----+----+----+----+----+----+
	  | B* | R* | B* | R* | B* | R* |
	  +----+----+----+----+----+----+
	  | R* | B* | R* | B* | R* | B* |
	  +----+----+----+----+----+----+
	  | B* | R* | B* | R* | B* | R* |
	  +----+----+----+----+----+----+
	*/
	shift(1);

        for(i = 1; i < m + 1; ++i)
        {
            for(j = jsta + (even(i) ^ 1); j <= jend; j += 2)
            {
 		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		diff = fabs(temp - uold[i][j]);
		if (diff > maxdiff) { maxdiff = diff; }
		uold[i][j] += omega * temp;
            }
        }

        ++iters;

	// Always check for symbolic difference.
	MPI_Allreduce(&maxdiff, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	maxdiff = diff;
    }
    while(maxdiff > stopdiff);
	    
    return iters;
}

// Divide portion of matrix according to rank and number of processors.
void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict jsta, int * restrict jend)
{
    const int iwork = ((n2 - n1) / nprocs) + 1;
    
    *jsta = min(irank * iwork + n1, n2 + 1);
    *jend = min(*jsta + iwork - 1, n2);
}

// Core MPI routine for red black SOR.
void shift(const int iflg)
{
    // Send and receive buffers.
    double bufs1[M + 2] = {0};
    double bufs2[M + 2] = {0};
    
    double bufr1[M + 2] = {0};
    double bufr2[M + 2] = {0};

    // Send and receive ranges.
    const int is1 = ((jsta + iflg) % 2) + 1;
    const int is2 = ((jend + iflg) % 2) + 1;

    const int ir1 = fabs(3 - is1);
    const int ir2 = fabs(3 - is2);

    int i, icnt1=0, icnt2=0;

    // Fill buffer to send to iprev
    if(rank != 0)
    {
	icnt1 = 0;

	for(i = is1; i < M + 2; i += 2)
	{
	    bufs1[icnt1] = uold[i][jsta];
	    ++icnt1;
	}
    }

    // Fill buffer to send to inext
    if (rank != numtasks - 1)
    {
	icnt2 = 0;

#ifdef DEBUG
	printf("rank = %d, ied=%d\n", rank, jend);
#endif
	
	for(i = is2; i < M + 2; i += 2)
	{
	    bufs2[icnt2] = uold[i][jend];

#ifdef DEBUG
	    printf("s%f ", bufs2[icnt2]);
#endif
	    ++icnt2;
	}
    }

    // Send data to iprev and inext
    MPI_Isend(bufs1, icnt1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqs1);
    MPI_Isend(bufs2, icnt2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqs2);
    
    // Receive data from iprev and inext
    MPI_Irecv(bufr1, M + 2, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqr1);
    MPI_Irecv(bufr2, M + 2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqr2);
    
    // Wait until everyone is finished
    MPI_Wait(&ireqs1, &istatus);
    MPI_Wait(&ireqs2, &istatus);
    
    MPI_Wait(&ireqr1, &istatus);
    MPI_Wait(&ireqr2, &istatus);

    int icnt;
   
    // Receive buffers from iprev
    if(rank != 0)
    {
	icnt = 0;

#ifdef DEBUF
	printf("rank = %d, jsta=%d\n", rank, jsta);
#endif

	for(i = ir1; i < M + 2; i += 2)
	{
	    uold[i][jsta - 1] = bufr1[icnt];
	    
#ifdef DEBUG
	    printf("r%f ", bufr1[icnt]);
#endif
	    ++icnt;
	}
    }

    // Receive buffers from inext
    if(rank != numtasks - 1)
    {
	icnt = 0;

	for(i = ir2; i < M + 2; i += 2)
	{
	    uold[i][jend + 1] = bufr2[icnt];
	    ++icnt;
	}
    }
}

inline int min(const int a, const int b)
{
    if(a > b)
    {
	return b;
    }

    return a;
}
