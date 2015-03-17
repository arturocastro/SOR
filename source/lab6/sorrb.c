#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
 
#define M 100

double compute_error(double solution[][M + 2], double u[][M + 2], const int m);
int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, const int m);

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
    }

    for(i = 0; i < m + 2; ++i)
    {
	for(j = 0; j < m + 1; ++j)
        {
            uold[i][j] = j * h * uold[i][M + 1];
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

int sor(double unew[][M + 2], double uold[][M + 2], double solution[][M + 2], const double omega, const double tol, 
const int m)
{
    int i, j;

    // Copy boudary conditions.
    for(i = 0; i < m + 2; ++i)
    {
        unew[i][m + 1] = uold[i][m + 1];
        unew[m + 1][i] = uold[m + 1][i];
        unew[i][0] = uold[i][0];
        unew[0][i] = uold[0][i];
    }

    int iters = 0;
    double error = compute_error(solution, uold, m);

    // Do SOR until 'tol' is satisfied
    while(error > tol)
    {
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
        for(i = 1; i < m + 1; ++i)
        {
            for(j = (i % 2) + 1; j < m + 1; j += 2)
            { 
                unew[i][j] = uold[i][j] + 0.25 * omega *
                                 (unew[i - 1][j] + unew[i][j - 1] +
                                 uold[i + 1][j] + uold[i][j + 1] -
                                 4.0 * uold[i][j]);
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
        for(i = 1; i < m + 1; ++i)
        {
            for(j = ((i + 1) % 2) + 1; j < m + 1; j += 2)
            { 
                unew[i][j] = uold[i][j] + 0.25 * omega *
                                 (unew[i - 1][j] + unew[i][j - 1] +
                                 uold[i + 1][j] + uold[i][j + 1] -
                                 4.0 * uold[i][j]);
            }
        }

	// Copy new to old.
        for(i = 1; i < m + 1; ++i)
        {
	    memcpy(&uold[i][1], &unew[i][1], m * sizeof(double));
        }

        ++iters;

	// Check every 20 iterations.
        if(iters % 20 == 0)
        {
            error = compute_error(solution, uold, m);
        }
    }

    return iters;
}
