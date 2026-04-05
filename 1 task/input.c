#include "input.h"

int input_mpi(int s, char* filename, int n, double* local_a, double* buffer, int rank, int size) {
	if (s == 0) return finput_mpi(filename, n, local_a, buffer, rank, size);
	else return sinput_mpi(s, n, local_a, rank, size);
}

int finput_mpi(char* filename, int n, double* local_a, double* buffer, int rank, int size) {
	int i = 0, j = 0, flag = 0; 
	FILE* in = NULL; 
	double curr = 0;

	if (rank == 0) {
		in = fopen(filename, "r");
		if (in == NULL) flag = 1;
	}
	MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (flag) return -1;

	for (i = 0; i < n; ++i) {
		if (rank == 0) {
			for (j = 0; j < n; ++j) {
				if (fscanf(in, "%lf", &curr) < 1) {
					flag = 1;
					break;
				}
				else buffer[j] = curr;
			}
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (flag) {
			if (rank == 0) fclose(in);
			return -2;
		}
	
		if (rank == 0) {
			if (i%size == 0) {
				for (j = 0; j < n; ++j) local_a[(i/size)*n+j] = buffer[j];
			} else {
				MPI_Send(buffer, n, MPI_DOUBLE, i%size, 0, MPI_COMM_WORLD);
			}
		} else if (rank == i%size) {
			//i = k*size + rank
			MPI_Recv(local_a+((i-rank)/size)*n, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	if (rank == 0) fclose(in);
	return 0;
}

int sinput_mpi(int s, int n, double* local_a, int rank, int size) {
	int i = 0, j = 0, local_n = 0;
	
	if (s < 1 || s > 4) {return -3;}
	
	local_n = n/size;
	if (rank < n%size) ++local_n;
	
	for (i = 0; i < local_n; ++i) {
		for (j = 0; j < n; ++j) {
			local_a[i*n+j] = formula_mpi(s, n, i*size+rank, j);
		}
	}

	return 0;
}

//formuls changed because (i, j) starts from (0, 0)
double formula_mpi(int s, int n, int i, int j) {
	switch (s) {
		case 1:
			return n- (i > j ? i : j);
		case 2:
			return (i > j ? i : j)+1;
		case 3:
			return (i > j ? i : j) - (i < j ? i : j);
		case 4:
			return 1./(i+j+1);
	}
	return 0;
}
