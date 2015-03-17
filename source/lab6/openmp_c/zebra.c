#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

double compute_error(double solution[][M + 2], double u[][M + 2], const int m);
int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m);

void printmat(double uold[][M + 2])
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

int main(int argc, char * argv[])
{
    const int m = M;

    printf("%d\n", m);

    double unew[M + 2][M + 2] = {{ 0 }};
    double solution[M + 2][M + 2] = {{ 0 }};

    double uold[M + 2][M + 2] = {{ 0 }};

    int i, j;

    const double begin = omp_get_wtime();

    const double pi = 4.0 * atan(1.0);

    const double h = pi / (m + 1);

    for(i = 0; i < m + 2; ++i)
    {
	uold[i][M + 1] = sin(i * h);
//	printf("%f\n", uold[i][M + 1]);
    }

//    exit(1);

    for(i = 0; i < m + 2; ++i)
    {
	for(j = 0; j < m + 1; ++j)
        {
            uold[i][j] = j * h * uold[i][M + 1];
//	    printf("%.40f\n", uold[i][j]);
        }
    }

//    exit(1);

    for(i = 0; i < m + 2; ++i)
    {
        for(j = 0; j < m + 2; ++j)
        {
            solution[i][j] = sinh(j * h) * sin(i * h) / sinh(pi);
	    // printf("%.40f\n", solution[i][j]);
        }
    }

    const double omega = 2.0 / ( 1.0 + sin(pi / (m + 1)) );
    const double tol = 0.001;

    const int iters = sor(unew, uold, solution, omega, tol, m);

    const double end = omp_get_wtime();

    printf(" \n");
    printf(" Omega = %f\n", omega);
    printf(" It took %d iterations.\n", iters);

    printf("Total time = %f\n\n\n", end - begin);

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

int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m)
{
    int i, j;
   
    /*for(i = 0; i < m + 2; ++i)
    {
        unew[i][m + 1] = uold[i][m + 1];
        unew[m + 1][i] = uold[m + 1][i];
        unew[i][0] = uold[i][0];
        unew[0][i] = uold[0][i];
    }*/

    double temp;
    int iters = 0;
    double error = compute_error(solution, uold, m);

    while(error > tol)
    {
#pragma omp parallel for private(i, j, temp)
        for(i = 1; i < m + 1; ++i)
        {
            for(j = (i % 2) + 1; j < m + 1; j += 2)
            {
//		printf("i=%d, j=%d\n", i, j);
//		printf("uold[i + 1][j]=%f\n"
		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		uold[i][j] += omega * temp;
		    //unew[i][j] = uold[i][j] + 0.25 * omega *
		    //           (unew[i - 1][j] + unew[i][j - 1] +
		    //           uold[i + 1][j] + uold[i][j + 1] -
		    //           4.0 * uold[i][j]);
		//printf("%f\n", unew[i][j]);
            }
        }

//	printmat(uold);
//	exit(1);

#pragma omp parallel for private(i, j, temp)
        for(i = 1; i < m + 1; ++i)
        {
            for(j = ((i + 1) % 2) + 1; j < m + 1; j += 2)
            {
		//printf("i=%d, j=%d\n", i, j);
		temp = 0.25 * (uold[i][j - 1] + uold[i - 1][j]
			       + uold[i + 1][j] + uold[i][j + 1]) - uold[i][j];
		uold[i][j] += omega * temp;
                //unew[i][j] = uold[i][j] + 0.25 * omega *
		//               (unew[i - 1][j] + unew[i][j - 1] +
		//               uold[i + 1][j] + uold[i][j + 1] -
		//               4.0 * uold[i][j]);
		//printf("%f\n", unew[i][j]);
            }
        }

	/*for(i = 1; i < m + 1; ++i)
	{
           for(j = 1; j < m + 1; ++j)
            { 
		printf("%f\n", unew[i][j]);
            }
        }*/

//	exit(1);

//        for(i = 1; i < m + 1; ++i)
//        {
//	    memcpy(&uold[i][1], &unew[i][1], m * sizeof(double));
//        }

        ++iters;

        if(iters % 20 == 0)
        {
            error = compute_error(solution, uold, m);
//	    printf("%f\n", error);
        }

//	if(iters == 42)
//	    exit(1);
    }

    return iters;
}
