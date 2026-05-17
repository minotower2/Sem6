#ifndef JORDAN_REVERSE_MPI
#define JORDAN_REVERSE_MPI

#include <stdlib.h>
#include <math.h>
#include <mpi.h>

typedef struct {
	double value;
	int index;
} MAX;

int reflect_reverse_mpi(int, double*, double*, int, int, double*, double*);

void r_mpi(int, double*, double*, double*, int, int, double*); 
void r1_r2_mpi(int, double*, double*, double*, double*, int, int, double*);

#endif
