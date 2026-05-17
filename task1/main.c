#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "input.h"
#include "output.h"
#include "reflect_reverse_MPI.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

	int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	int local_n = 0, n = 0, r = 0, s = 0, task = 0, flag = 0, i = 0, k = 0, j = 0;
	char *filename = NULL; 
	double *local_array = NULL, *local_rev = NULL, *buffer = NULL, *local_d = NULL; 
	double r1 = 0, r2 = 0, t1 = 0, t2 = 0;
	
	task = 24;
	//Initialization
	if (argc < 4 || argc > 5) {
		printf("Invalid number of arguments.\n");
		MPI_Finalize();
		return -1;
	}
	n = strtol(argv[1], NULL, 10);
	r = strtol(argv[2], NULL, 10);
	s = strtol(argv[3], NULL, 10);

	if (s == 0 && argc == 4) {
		printf("Invalid number of arguments.\n");
		MPI_Finalize();
		return -1;
	}

	local_n = n/size;
	if (rank < n%size) ++local_n;
	
	local_array = (double*)malloc((n/size+1)*n*sizeof(double));
	local_d = (double*)malloc((n/size+1)*sizeof(double));


	if (rank == 0 && s == 0) {
		filename  = argv[4]; 
	}
	buffer = (double*)malloc(2*n*sizeof(double));

	flag = input_mpi(s, filename, n, local_array, buffer, rank, size);
	switch (flag) {
		case -1:
			if(rank == 0) printf("Error during opening file %s\n", filename);
			free(local_array);
			free(buffer);
			MPI_Finalize();
			return -1;
		case -2:
			if (rank == 0) printf("Wrong data in file %s\n", filename);
			free(local_array);
			free(buffer);
			MPI_Finalize();
			return -1;
		case -3:
			if (rank == 0) printf("Wrong type of formule: %d\n", s);
			free(local_array);
			free(buffer);
			MPI_Finalize();
			return -1;
	}

	output_file_mpi("out_data.txt", r, r, n, local_array, rank, size);	
	
	local_rev = (double*)malloc((n/size + 1)*n*sizeof(double));

	for (i = 0; i < local_n; ++i) {
		for (j = 0; j < n; ++j) {
			if (j == i*size+rank)
				local_rev[i*n+j] = 1;
			else 
				local_rev[i*n+j] = 0;
		}
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	flag = reflect_reverse_mpi(n, local_array, local_rev, rank, size, local_d, buffer);
	
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime()-t1;

	if (flag) {	
		r1 = -1; r2 = -1; t2 = 0;
	} else {
		output_file_mpi("out_rev.txt", r, r, n, local_rev, rank, size);	

		flag = input_mpi(s, filename, n, local_array, buffer, rank, size);

		MPI_Barrier(MPI_COMM_WORLD);
		t2 = MPI_Wtime();

		r1_r2_mpi(n, local_array, local_rev, &r1, &r2, rank, size, buffer);
	
		MPI_Barrier(MPI_COMM_WORLD);
		t2 = MPI_Wtime()-t2;
	}
	
	if (rank == 0) {
		printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d\n", argv[0], task, r1, r2, t1, t2, s, n);
	}
	
	free(local_array);
	free(local_rev);
	free(buffer);

	MPI_Finalize();
	return 0;
}
