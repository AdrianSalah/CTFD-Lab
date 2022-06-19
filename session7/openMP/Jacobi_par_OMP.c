/******************************************************************************
* FILE: Jacobi_par_OMP.c
* DESCRIPTION:  Jacobi algorithm parallelized using openMP
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*** Define matrix size*/

#define N 1000           	/* Size of matrix*/
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
int	iter, i, j, chunk;

omp_set_num_threads(10); /*Set number of threads */

epsilon = epsilon*epsilon;

/*** Initialisation */
gettimeofday(&start, NULL);
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

/*** Computation loop ***/
/* Parallel section */
chunk = omp_get_max_threads();
printf("Nb processors : %i \n",omp_get_num_procs());
printf("Nb threads : %i\n",chunk);
chunk = (int)((N)/chunk/2);
printf("Chunksize: %i\n",chunk);

/*
// Get environment information
int nthreads, tid, procs, maxt, inpar, dynamic, nested;
procs = omp_get_num_procs();
nthreads = omp_get_num_threads();
maxt = omp_get_max_threads();
inpar = omp_in_parallel();
dynamic = omp_get_dynamic();
nested = omp_get_nested();

// Print environment information 
printf("Nb processors = %d\n", procs);
printf("Nb of threads = %d\n", nthreads);
printf("Max threads = %d\n", maxt);
printf("In parallel? = %d\n", inpar);
printf("Dynamic threads enabled? = %d\n", dynamic);
printf("Nested parallelism supported? = %d\n", nested);  
*/

gettimeofday(&start, NULL);
for(iter=0; conv > epsilon; iter++){
	conv = 0;
	// Parallel part
	#pragma omp parallel default(none) shared(A,b,x1,x2,conv,chunk) private(i,j)
	{
		#pragma omp for //schedule (dynamic,chunk)
		for(i = 0; i<N; i++){
			x2[i] = 0;
			for(j=0; j<N; j++){
				x2[i] += A[i][j]*x1[j];
			}
			x2[i]= (b[i]-x2[i])/A[i][i] + x1[i];
		}
		#pragma omp for //schedule (dynamic,chunk)
		for (i=0; i<N; i++){
			#pragma omp reduction(+:conv)
			{
			conv += (x2[i] - x1[i])*(x2[i] - x1[i]);
			}
			x1[i] = x2[i];
		}
	}
}
gettimeofday(&end, NULL);

/*** Result print */

cpu_time_used = (double)((end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec))/1000000.0;

printf("Conv^2: %E ;\nIter: %d ;\nTime: %f ;\nInit Time: %f \n",conv,iter,cpu_time_used,cpu_time_init);
printf("%f \t %f \t %f \t %f \n",x2[0],x2[1],x2[2],x2[3]);

return 0;
}
