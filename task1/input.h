#ifndef INPUT_MPI
#define INPUT_MPI

#include <stdio.h>
#include <mpi.h>

int input_mpi(int, char*, int, double*, double*, int, int);

int finput_mpi(char*, int, double*, double*, int, int);

int sinput_mpi(int, int, double*, int, int);
double formula_mpi(int, int, int, int);

#endif
