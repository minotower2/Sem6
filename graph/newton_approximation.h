#ifndef NEWTON_APPROXIMATION_H
#define NEWTON_APPROXIMATION_H

void make_Lagrange_polynomial(int n, double* x, double* f, double* d, double* c);
double calculate_newton_approximation(double x_0, int n, double* x, double* c);
double calculate_newton_discrepancy(double x_0, func_t func, int n, double* x, double* c);

#endif