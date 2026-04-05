#include <stdio.h>
#include "output.h"
#include "reflect_reverse_MPI.h"

int reflect_reverse_mpi(int n, double *local_a, double *local_rev, int rank, int size, double *local_d, double* buffer) {
	int local_n = 0, k = 0, start = 0, local_n_max = 0, j = 0, i = 0, t = 0;
	double local_max = 0, tmp = 0, buff = 0, local_s = 0, s = 0, norm_a1 = 0, prod = 0;
	
	local_n = n/size;
	if (rank < n%size)
		++local_n;
	
	for (k = 0; k < n; ++k) {
		s = 0;
		local_s = 0;
		//i*size+rank >= k+1
		for (i = 0; i < local_n; ++i) {
			if (i*size+rank >= k+1)
				local_s += local_a[i*n+k]*local_a[i*n+k];
		}

		MPI_Reduce(&local_s, &s, 1, MPI_DOUBLE, MPI_SUM, k%size, MPI_COMM_WORLD);

		if (rank == k%size) {
			norm_a1 = sqrt(local_a[((k-rank)/size)*n+k]*local_a[((k-rank)/size)*n+k]+s);
			local_a[((k-rank)/size)*n+k] -= norm_a1;
			local_d[(k-rank)/size] = norm_a1;
			s = sqrt(local_a[((k-rank)/size)*n+k]*local_a[((k-rank)/size)*n+k]+s);
		}

		MPI_Bcast(&norm_a1, 1, MPI_DOUBLE, k%size, MPI_COMM_WORLD);

		if (norm_a1 < 1e-15)
			return -1;

		MPI_Bcast(&s, 1, MPI_DOUBLE, k%size, MPI_COMM_WORLD);

		if (s < 1e-15)
			continue;

		for (i = 0; i < local_n; ++i) {
			if (i*size+rank >= k)
				local_a[i*n+k] /= s;
		}

		for (i = k+1; i < n; ++i) {
			buff = 0;
			for (t = 0; t < local_n; ++t) {
				if (t*size+rank >= k)
					buff += local_a[t*n+k]*local_a[t*n+i];
			}

			MPI_Allreduce(&buff, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			prod *= 2;

			for (t = 0; t < local_n; ++t) {
				if (t*size+rank >= k)
					local_a[t*n+i] -= prod*local_a[t*n+k];
			}

		}

		for (i = 0; i < n; ++i) {
			buff = 0;
			for (t = 0; t < local_n; ++t) {
				if (t*size+rank >= k)
					buff += local_a[t*n+k]*local_rev[t*n+i];
			}

			MPI_Allreduce(&buff, &prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			prod *= 2;

			for (t = 0; t < local_n; ++t) {
				if (t*size+rank >= k)
					local_rev[t*n+i] -= prod*local_a[t*n+k];
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	for (k = n-1; k >= 0; --k) {
		if (rank == k%size) {
			for (j = 0; j < n; ++j) 
				buffer[j] = local_a[((k-rank)/size)*n+j];
		}

		MPI_Bcast(buffer, n, MPI_DOUBLE, k%size, MPI_COMM_WORLD);

		for (i = 0; i < n; ++i) {
			buffer[i+n] = 0;
			for (j = 0; j < local_n; ++j) {
				if (j*size+rank >= k+1)
					buffer[i+n] += buffer[j*size+rank]*local_rev[j*n+i];
			}
		}

		MPI_Reduce(buffer+n, buffer, n, MPI_DOUBLE, MPI_SUM, k%size, MPI_COMM_WORLD);
		if (rank == k%size) {
			for (i = 0; i < n; ++i)
				local_rev[(k-rank)/size*n+i] = (local_rev[(k-rank)/size*n+i] - buffer[i])/local_d[(k-rank)/size];
		}
	}

	return 0;
}

void r_mpi(int n, double* local_a, double* local_b, double* r, int rank, int size, double* buffer) {
	double norm_ab_e = 0, sum = 0, max = 0;
	int i = 0, k = 0, new_rank = 0, tmp_n = 0, local_n = 0, j = 0, t = 0;
	MPI_Datatype local_col_t;

	MPI_Type_vector(n/size+1, 1, n, MPI_DOUBLE, &local_col_t);
	MPI_Type_commit(&local_col_t);

	norm_ab_e = -1;

	local_n = n/size;
	if (rank < n%size) 
		++local_n;

	for (t = 0; t < n; ++t) {
		sum = 0;
		for (i = 0; i < local_n; ++i) {
			if (i*size+rank == t)
				buffer[i] = -1;
			else
				buffer[i] = 0;
		}

		for (k = 0; k < size; ++k) {
			new_rank = (rank+k)%size;
			tmp_n = n/size;
			if (new_rank < n%size) 
				++tmp_n;

			for (i = 0; i < local_n; ++i) {
				for (j = 0; j < tmp_n; ++j) {
					buffer[i] += local_a[i*n+(j*size+new_rank)]*local_b[j*n+t];
				}
			}	
		
			MPI_Sendrecv_replace(local_b+t, 1 , local_col_t, (rank - 1 + size) % size, 0, (rank+1)%size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		for (i = 0; i < local_n; ++i) {
			sum += fabs(buffer[i]);
		}
		MPI_Reduce(&sum, &max, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0 && norm_ab_e < max) {
			norm_ab_e = max;
		}
	}
	MPI_Type_free(&local_col_t);
	if (rank == 0)
		*r = norm_ab_e;
} 

void r1_r2_mpi(int n, double* local_a, double* local_rev, double* r1, double* r2, int rank, int size, double* buffer) {
	if (n > 11000) {
		*r1 = 0;
		*r2 = 0;
		return;
	}

	r_mpi(n, local_a, local_rev, r1, rank, size, buffer);
	r_mpi(n, local_rev, local_a, r2, rank, size, buffer);
}
