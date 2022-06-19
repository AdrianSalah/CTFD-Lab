/******************************************************************************
* FILE: Jacobi_par_OMP.c
* DESCRIPTION:  Jacobi algorithm parallelized using MPI
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*** Define matrix size*/

#define N 1000      /* Size of the matrix */
#define diag 1000	/* Value on the diagonal */
#define block 1		/* Only 1 block */

/*** Parameters declaration initialization */

double	x1[N],			/* Iterative vector */
	epsilon=1E-5,		/* Convergence criteria */
	conv=1;			/* Convergence values */
int	subN;			/* size of subblock */
int	maxIter = 20000;	/* Max number of iteration */
double 	cpu_time_used,		/* Computation time */
	cpu_time_init;		/* Initialisation time */
struct 	timeval start, end;

int main (int argc, char *argv[]) 
{
/*** Initialisation ***/
int iter, i, j, rank, nbSlave;
int stop =1; // Stop variable; value 0 sent to all slaves"

int	Amaster[N*N],	/* Declaration A */
	bmaster[N];     /* Declaration B */
double	x2master[N];
	
/*** Computation loop ***/
MPI_Status status;
MPI_Datatype Amatrix;
MPI_Datatype bvec;
MPI_Datatype x1vec;
MPI_Datatype x2vec;

MPI_Init(&argc, &argv); // Initialize the MPI environment
MPI_Comm_rank(MPI_COMM_WORLD,&rank); // get the thread nb in rank
MPI_Comm_size(MPI_COMM_WORLD,&nbSlave); // get the number of workers

subN = N/(nbSlave-1); // Each worker computes the blocksize

MPI_Type_vector(block,N*subN,N*subN,MPI_INT,&Amatrix); // array of size N*subN, to send/receive part of A

MPI_Type_vector(block,subN,subN,MPI_INT,&bvec); // array of size subN to send/receive part of b
MPI_Type_vector(block,N,N,MPI_DOUBLE,&x1vec); // array of size N to send/receive x1
MPI_Type_vector(block,subN,subN,MPI_DOUBLE,&x2vec); // array of size subN to send/receive the parts of x2

MPI_Type_commit(&Amatrix);
MPI_Type_commit(&bvec);
MPI_Type_commit(&x1vec);
MPI_Type_commit(&x2vec);

// Local slave variables only
double x2[subN]; // iterative parts of x2
int b[subN];
int A[subN*N]; // Subpart of A

if(rank==0) // master only
{
	printf("start calculation \n");

	epsilon = epsilon*epsilon;
	gettimeofday(&start, NULL);

	printf("Initialization\n");

	for(i=0; i < N*N; i++)
		Amaster[i] = 1;
		
	for (i=0; i < N; i++) {
		bmaster[i] = 1;
		x1[i] = 1;
		Amaster[i+i*N] = diag;
	}
	
	gettimeofday(&end, NULL);
	cpu_time_init = (double)((end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec))/1000000.0;

	printf("\t Initialization time : %f\n", cpu_time_init);
	
	gettimeofday(&start, NULL);
	// Send the parts of A and b to each slaves
	for ( i=1; i < nbSlave;i++) {
		MPI_Send(&Amaster[subN*(i-1)*N],1,Amatrix,i,i,MPI_COMM_WORLD);
		MPI_Send(&bmaster[(i-1)*subN],1,bvec,i,i,MPI_COMM_WORLD);
	}
	
	for (iter = 1; iter < maxIter && epsilon < conv ; iter++) {
		conv = 0;
		MPI_Bcast(&x1,1,x1vec,0,MPI_COMM_WORLD); // send x1 to everyone
			
		for (i=1 ; i < nbSlave;i++){
			MPI_Recv(&x2master[(i-1)*subN],subN,MPI_DOUBLE,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			//MPI_Recv(&convIter,1,MPI_DOUBLE,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			//conv += convIter;
		}
		
		// only the master computes the convergence (otherwise too much send/receive)
		for (i=0; i < N	; i++) {
			conv += (x2master[i] - x1[i])*(x2master[i] - x1[i]);
			x1[i] = x2master[i];
		}
		
		if (conv < epsilon || iter == maxIter-1)
			stop = 0;
		
		// send of stop signal to everyone
		MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	gettimeofday(&end, NULL);
	cpu_time_used = (double)((end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec))/1000000.0;
	printf("Conv^2: %E ;\nIter: %d ;\nTime: %f ;\nInit Time: %f \n",conv,iter,cpu_time_used,cpu_time_init);
	printf("%f \t %f \t %f \t %f \n",x1[0],x1[1],x1[2],x1[3]);
		
} else // slave work
{
	MPI_Recv(&A[0],1,Amatrix,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	MPI_Recv(&b[0],1,bvec,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	while(stop)
	{
		MPI_Bcast(&x1,1,x1vec,0,MPI_COMM_WORLD);
		
		conv = 0;
		for (i=0; i < subN; i++) {
			x2[i] = 0;
			for(j=0; j<N; j++){
				x2[i] += A[i*N+j]*x1[j];
			}
			x2[i] = x1[i+(rank-1)*subN] + (b[i] - x2[i])/A[(i)*N+(rank-1)*subN+i];
			//conv += (x2[i] - x1[i+(rank-1)*subN])*(x2[i] - x1[i+(rank-1)*subN]);
		}
		MPI_Send(&x2[0],1,x2vec,0,rank,MPI_COMM_WORLD);
		//MPI_Send(&conv,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		
		MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
}
MPI_Finalize();
return 0;
}

