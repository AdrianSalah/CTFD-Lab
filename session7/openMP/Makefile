CC=gcc
CFLAGS=-fopenmp

% : %c
	$(CC) $(CFLAGS) -o % $<

Jacobi_par_OMP : Jacobi_par_OMP.c
	$(CC) $(CFLAGS) -o Jacobi_par_OMP Jacobi_par_OMP.c
