#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
 
#define M 20

int numtasks, rank, inext, iprev, ista, iend;

double unew[M + 2][M + 2] = {{ 0 }};
double solution[M + 2][M + 2] = {{ 0 }};

double uold[M + 2][M + 2] = {{ 0 }};

MPI_Request ireqs1;
MPI_Request ireqs2;

MPI_Request ireqr1;
MPI_Request ireqr2;

MPI_Status istatus;

double compute_error(double solution[][M + 2], double u[][M + 2], const int m);
int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m);
void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict ista, int * restrict iend);
void shift(const int iflg);
inline int min(const int a, const int b);

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
    for (i = ista; i <= iend; ++i)
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

    for(i = 0; i < m + 2; ++i)
    {
        for(j = 0; j < m + 2; ++j)
        {
            solution[i][j] = sinh(j * h) * sin(i * h) / sinh(pi);
        }
    }

    para_range(1, m + 1, numtasks, rank, &ista, &iend);

    printf("numtasks = %d\nrank %d: from %d to %d\n", numtasks, rank, ista, iend);

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

    const int iters = sor(unew, uold, solution, omega, tol, m);

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

double compute_error(double solution[][M + 2], double u[][M + 2], const int m)
{
    double error = 0.0;
    int i, j;

    for(i = ista; i <= iend; ++i)
    {
        for(j = 1; j < m + 1; ++j)
        {
	    const double abs_diff = fabs(solution[i][j] - u[i][j]);

            if(error < abs_diff)
	    {
                error = abs_diff;
	    }
        }
    }

    return error;
}

int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m)
{
    int i, j;

    int iters = 0;
    double error = compute_error(solution, uold, m);
    double error2 = error;
    double temp;

    while(error2 > tol)
    {
#ifdef DEBUG
	if(rank == 1)
	{
	    printf("ista = %d, iend = %d\n", ista, iend);

	    printmat();	    
	}
#endif

	shift(0);

#ifdef DEBUG
	if(rank == 0){puts("BOOM!");printmat();}
#endif

        for(i = ista; i < iend; ++i)
        {
            for(j = (i % 2) + 1; j < m + 1; j += 2)
            {
 		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		uold[i][j] += omega * temp;
            }
        }

	/*if (rank == 1)
	{
	    sleep(1);
	}
	printmatrank();
	sleep(5);
	exit(1);*/

//	if(rank == 0){puts("BOOoM!");printmat();}
	shift(1);
//	if(rank == 1){puts("BOOooM!");printmat();}

        for(i = ista; i < iend; ++i)
        {
            for(j = ((i + 1) % 2) + 1; j < m + 1; j += 2)
            {
		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		uold[i][j] += omega * temp;
            }
        }

        ++iters;

	//exit(1);

        if(iters % 20 == 0)
        {
            error = compute_error(solution, uold, m);
	    MPI_Allreduce(&error, &error2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    printf("%d=%f\n", rank, error2);
        }
    }

    return iters;
}

void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict ista, int * restrict iend)
{
    const int iwork = ((n2 - n1) / nprocs) + 1;
    
    *ista = min(irank * iwork + n1, n2 + 1);
    *iend = min(*ista + iwork - 1, n2);
}

void shift(const int iflg)
{
    double bufs1[M + 2] = {0};
    double bufs2[M + 2] = {0};
    
    double bufr1[M + 2] = {0};
    double bufr2[M + 2] = {0};

    const int is1 = ((ista + iflg) % 2) + 1;
    const int is2 = ((iend + iflg) % 2) + 1;

    const int ir1 = fabs(3 - is1);
    const int ir2 = fabs(3 - is2);

    int i, icnt1=0, icnt2=0;

    if(rank != 0)
    {
	icnt1 = 0;

	for(i = is1; i < M + 1; i += 2)
	{
	    bufs1[icnt1] = uold[ista][i];
	    ++icnt1;
	}
    }

    if (rank != numtasks - 1)
    {
	icnt2 = 0;

#ifdef DEBUG
	printf("rank = %d, ied=%d\n", rank, iend);
#endif
	
	for(i = is2; i < M + 1; i += 2)
	{
	    bufs2[icnt2] = uold[iend][i];

#ifdef DEBUG
	    printf("s%f ", bufs2[icnt2]);
#endif
	    ++icnt2;
	}
    }

    //if(iflg == 0)
	//printf("iflg=%d\nis1 = %d, is2 = %d\nir1 = %d, ir2=%d\n", iflg, is1, is2, ir1, ir2);

    MPI_Isend(bufs1, icnt1, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqs1);
    MPI_Isend(bufs2, icnt2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqs2);
    
    MPI_Irecv(bufr1, M + 2, MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &ireqr1);
    MPI_Irecv(bufr2, M + 2, MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &ireqr2);
    
    MPI_Wait(&ireqs1, &istatus);
    MPI_Wait(&ireqs2, &istatus);
    
    MPI_Wait(&ireqr1, &istatus);
    MPI_Wait(&ireqr2, &istatus);

    int icnt;
   
    if(rank != 0)
    {
	icnt = 0;

#ifdef DEBUF
	printf("rank = %d, ista=%d\n", rank, ista);
#endif

	for(i = ir1; i < M + 1; i += 2)
	{
	    uold[ista - 1][i] = bufr1[icnt];
	    
#ifdef DEBUG
	    printf("r%f ", bufr1[icnt]);
#endif
	    ++icnt;
	}
    }

    if(rank != numtasks - 1)
    {
	icnt = 0;

	for(i = ir2; i < M + 1; i += 2)
	{
	    uold[iend + 1][i] = bufr2[icnt];
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
