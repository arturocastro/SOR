#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
 
#define M 500

double compute_error(double solution[][M + 2], double u[][M + 2], const int m);
int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m, const double h2);

int main(void)
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

    struct timeval tic, toc;

    double unew[M + 2][M + 2] = {{ 0 }};
    double solution[M + 2][M + 2] = {{ 0 }};

    double uold[M + 2][M + 2] = {{ 0 }};

    int i, j;

    gettimeofday(&tic, NULL);
    // const double begin = omp_get_wtime();

    const double pi = 4.0 * atan(1.0);

    const double h = pi / (m + 1);
    //const double h2 = (1 / (m + 1)) * (1 / (m + 1));

    for(i = 0; i < m + 2; ++i)
    {
       uold[i][M + 1] = sin(i * h);
    }


    for(i = 0; i < m + 2; ++i)
    {
	      for(j = 0; j < m + 1; ++j)
        {
            uold[i][j] = j * h * uold[i][M + 1];
//	    printf("%.40f\n", uold[i][j]);
        }
    }


    for(i = 0; i < m + 2; ++i)
    {
        for(j = 0; j < m + 2; ++j)
        {
            solution[i][j] = sinh(j * h) * sin(i * h) / sinh(pi);
//	    printf("%.40f\n", solution[i][j]);
        }
    }

    const double omega = 2.0 / ( 1.0 + sin(pi / (m + 1)) );
    const double tol = 0.001;

    const int iters = sor(unew, uold, solution, omega, tol, m, h2);

    gettimeofday(&toc, NULL);
    //const double end = omp_get_wtime();

    printf(" \n");
    printf(" Omega = %f\n", omega);
    printf(" It took %d iterations.\n", iters);

    printf("Total time = %f\n\n\n", (double) (toc.tv_usec - tic.tv_usec) / 1000000 +
	                            (double) (toc.tv_sec - tic.tv_sec));

    //printf("Total time = %f\n\n\n", end - begin);

    return 0;
}

double compute_error(double solution[][M + 2], double u[][M + 2], const int m)
{
    double error = 0.0;
    int i, j;

    for(i = 1; i < m + 1; ++i)
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

int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, 
	const int m, const double h2)
{
    int i, j;
   
    for(i = 0; i < m + 2; ++i)
    {
        unew[i][m + 1] = uold[i][m + 1];
        unew[m + 1][i] = uold[m + 1][i];
        unew[i][0] = uold[i][0];
        unew[0][i] = uold[0][i];
    }

    int iters = 0;
    double error = compute_error(solution, uold, m);

    while(error > tol)
    {
          for(i = 1; i < m + 1; ++i)
          {
              for(j = 1; j < m + 1; ++j)
              { 
                  unew[i][j] = (uold[i+1][j] + uold[i-1][j] +
				uold[i][j+1] + uold[i][j-1] - 
				( uold[i][j])) * 0.25;
              }
          }
        

        for(i = 1; i < m + 1; ++i)
        {
	        memcpy(&uold[i][1], &unew[i][1], m * sizeof(double));

        }

        ++iters;

        if(iters % 20 == 0)
        {
            error = compute_error(solution, uold, m);
        }
    }

    return iters;
}
