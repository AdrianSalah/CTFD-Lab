CC=mpicc
CFLAGS=-g

% : %c
	$(CC) $(CFLAGS) -o % $<

Jacobi_par_MPI : Jacobi_par_MPI.c
	$(CC) -DVERIFY -o Jacobi_par_MPI Jacobi_par_MPI.c
