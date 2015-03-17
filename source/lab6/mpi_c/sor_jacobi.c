#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
 
#define M 100

double compute_error(double solution[][M + 2], double u[][M + 2], const int m);
int sor(double uold[][M + 2],double unew[][M+2], double solution[][M+2],const double tol, const int m, const double h2);
void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict jsta, int * restrict jend);
void shift();
inline int min(const int a, const int b);

int numtasks, rank, inext, iprev, jsta, jend;

double uold[M + 2][M + 2] = {{ 0 }};
double unew[M + 2][M + 2] = {{ 0 }};

MPI_Request ireqs1;
MPI_Request ireqs2;

MPI_Request ireqr1;
MPI_Request ireqr2;

MPI_Status istatus;

void printmat(double mat[M+2][M+2])
{
    int i, j;
    for(i = 0; i < M + 2; ++i)
    {
	for(j = 0; j < M + 2; ++j)
	{
	    printf("%.2f ", mat[i][j]);
	}
	puts("\n");
    }

}
int main(int argc, char * argv[])
{

/*  
    Solution of Laplace's Equation.
    ==============================
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

*/
    const int m = M;


    int ierr = MPI_Init(&argc, &argv);

    if(ierr != MPI_SUCCESS)
    {
	perror("MPI init failed. Terminating T800.\n");
	exit(1);
    }

    double solution[M + 2][M + 2] = {{ 0 }};

    int i, j;

    const double tol = 0.001;
    
    const double begin = MPI_Wtime();

    const double pi = 4.0 * atan(1.0);

    const double h = pi / (m + 1);

    const double h2 = (1 / (m + 1)) * (1 / (m + 1));

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    //initialization (all the processes)
    for(i = 0; i < m + 2; ++i)
    {
       uold[i][M + 1] = sin(i * h);
    }

    for(i = 0; i < m + 2; ++i)
    {
	for(j = 0; j < m + 1; ++j)
        {
            uold[i][j] = j * h * uold[i][M + 1];
        }
    }


    for(i = 0; i < m + 2; ++i)
    {
        for(j = 0; j < m + 2; ++j)
        {
            solution[i][j] = sinh(j * h) * sin(i * h) / sinh(pi);
        }
    }
	
    //assign columns indexes to processes
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
    
    

    const int iters = sor(uold, unew, solution, tol, m, h2);

    const double end = MPI_Wtime();
    MPI_Finalize();

    if(rank == 0)
    {
     printf(" \n");
     printf(" Omega = %f\n", omega);
     printf(" It took %d iterations.\n", iters);

     printf("Total time = %f\n\n\n", end - begin);

     //printf("\n\nOld\n");

     //printmat(uold);

     //printf("\n\nSolution\n");

     //printmat(solution);
    }
    return 0;
}

double compute_error(double solution[][M + 2], double u[][M + 2], const int m)
{
    double error = 0.0;
    int i, j;

    for(i = 1; i < m + 1; ++i)
    {
        for(j = jsta; j <= jend; ++j)
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

int sor(double uold[][M + 2],double unew[][M+2], double solution[][M+2],const double tol, const int m, const double h2)
{
    int i, j;

    int iters = 0;
    //double error = compute_error(solution, uold, m);

    //initialize the unew matrix 
    for(i = 0; i < m + 2; ++i)
    {
        unew[i][m + 1] = uold[i][m + 1];
        unew[m + 1][i] = uold[m + 1][i];
        unew[i][0] = uold[i][0];
        unew[0][i] = uold[0][i];
    }

    //compute the error
    double error = compute_error(solution, uold, m);
    double global_error = error;
    do
    {
	

	//distribute and get the ghost cols
	shift();

        for(i = 1; i < m + 1; ++i)
        {
	    for(j = jsta; j <= jend; ++j)
            { 
		unew[i][j] = (uold[i+1][j] + uold[i-1][j] +
			      uold[i][j+1] + uold[i][j-1]- (h2*uold[i][j])  ) * 0.25;
				//- (h2*uold[i][j])F
            }
        }
        
	//copy the new into the old
        for(i = 1; i < m + 1; ++i)
        {
	    memcpy(&uold[i][1], &unew[i][1], m * sizeof(double));   
	}
	

	//error = compute_error(solution, uold, m);
	//if(iters % 20 == 0)
        //{
            error = compute_error(solution, uold, m);
	    //printf("process %d error %.10f global error %.10f \n", rank, error, global_error);
	    
        //}
	//printf("error   %.10f ", error);

        ++iters;
	//spread the local error maxdiff and get the others error in diff
	MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	error = global_error;
		

        
    }while(global_error > tol);

    return iters;
}

void para_range(const int n1, const int n2, const int nprocs, const int irank, int * restrict jsta, int * restrict jend)
{
    const int iwork = ((n2 - n1) / nprocs) + 1;
    
    *jsta = min(irank * iwork + n1, n2 + 1);
    *jend = min(*jsta + iwork - 1, n2);
}

void shift()
{
    // Send buffers.
    double bufs1[M + 2] = {0};
    double bufs2[M + 2] = {0};
    
    //receive buffer (to receive from the previous process the previous colums)
    double bufr1[M + 2] = {0};
    double bufr2[M + 2] = {0};


    int i, icnt1=0, icnt2=0;

    // Fill buffer to send to iprev
    if(rank != 0)
    {
		icnt1 = 0;

		for(i = 1; i < M + 2; i++)
		{
			bufs1[icnt1] = uold[i][jsta];
			++icnt1;
		}
    }

    // Fill buffer to send to inext
    if (rank != numtasks - 1)
    {
		icnt2 = 0;

	
		for(i = 1; i < M + 2; i++)
		{
			bufs2[icnt2] = uold[i][jend];

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


		//take jsta-1 the previous column from the previous process --OK
		for(i = 1; i < M + 2; i++)
		{
			uold[i][jsta - 1] = bufr1[icnt];

			++icnt;
		}
    }

    // Receive buffers from inext
    if(rank != numtasks - 1)
    {
		icnt = 0;
			//take jsta+1 the next column from the previous process --OK
			for(i = 1; i < M + 2; i++)
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

