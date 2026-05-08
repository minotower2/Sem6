#ifndef SPLINE_APPROXIMATION_H
#define SPLINE_APPROXIMATION_H

void make_spline(int n, double* x, double* f, double* c);
double calculate_spline_approximation(double x_0, int n, double* x, double* c);
double calculate_spline_discrepancy(double x_0, func_t func, int n, double* x, double* c);

#endif