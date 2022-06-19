/******************************************************************************
* FILE: Jacobi_seq.c
* DESCRIPTION:  Sequential version of Jacobi algorithm
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*** Define matrix size*/

#define N 1000            	/* Size of matrix */
#define diag 1000		/* Value on the diagonal */

/*** Parameters declaration initialization */

double	A[N][N],           	/* Declaration A */
	b[N],             	/* Declaration B */
	x1[N],x2[N],        	/* Iterative vector */
	epsilon=1E-4,		/* Convergence criteria */
	conv=1;			/* Convergence iteration */
int	maxIter= 100000;	/* Max number of iteration */
double 	cpu_time_used,		/* Computation time */
	cpu_time_init;
struct 	timeval start, end;

int main (int argc, char *argv[]) 
{
int	iter, i, j;

printf("Start calculation\n");

epsilon = epsilon*epsilon;
gettimeofday(&start, NULL);

printf("Initialization\n");

for(i=0; i<N; i++){
	for(j=0; j<N; j++){
		A[i][j] = 1;
	}
}
for(i=0; i<N; i++){
	A[i][i] = diag;
	b[i] = 1;
	x1[i] = 1;
}
gettimeofday(&end, NULL);
cpu_time_init = (double)((end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec))/1000000.0;

printf("\t Initialization time : %f\n", cpu_time_init);
printf("Start of calculation loop\n");

/*** Computation loop ***/
gettimeofday(&start, NULL);
for(iter=0; conv > epsilon; iter++){ // && iter < maxIter
	conv = 0;
	//printf("Iteration : %i\n",iter);
	for(i = 0; i<N; i++){
		x2[i] = 0;
		for(j=0; j<N; j++){
			x2[i] += A[i][j]*x1[j];
		}
		x2[i] = x1[i] + (b[i] - x2[i])/A[i][i];
		// conv += (x2[i] - x1[i])^2;
		conv += (x2[i] - x1[i])*(x2[i] - x1[i]);
	}
	for (i=0; i<N; i++)
		x1[i] = x2[i];
}
gettimeofday(&end, NULL);

/***  Show some results */

printf("End of computation \n");

cpu_time_used = (double)((end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec))/1000000.0;

printf("\tComputing time : %f\n", cpu_time_used);
printf("x[0] %f, x[1] %f, x[2] %f",x2[0],x2[1],x2[2]);
printf("\tConvergence^2 : %f\n\tIterations : %d\n",conv,iter);

return 0;
}
