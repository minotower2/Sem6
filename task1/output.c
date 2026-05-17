#include "output.h"

void output_file_mpi(char* filename, int r, int l, int n, double* local_a, int rank, int size) {
	int i = 0, j = 0;
	FILE *out = NULL;

	if (rank == 0) {
		out = fopen(filename, "w");
		fclose(out);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < (n < l ? n : l); ++i) {
		if (i%size == rank) {
			out = fopen(filename, "a");
			fprintf(out, "rank = %d       ", rank);
			for (j = 0 ; j < (r < n ? r : n); ++j) {
				fprintf(out, " %10.3e", local_a[(i-rank)/size*n+j]);
			}
			fprintf(out, "\n");
			fclose(out);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
