#ifndef PARABOLIC_APPROXIMATION_H
#define PARABOLIC_APPROXIMATION_H

#include "functions.h" // Для func_t

bool solve_tridiagonal(int n, const double* a, const double* b, const double* c,
                              const double* d, double* x);


int make_parabolic_spline_coefficients(int n,
                                       const double* x_nodes,
                                       const double* f_values,
                                       const double* xi_nodes,
                                       double* coeffs);


double calculate_parabolic_spline_approximation(double x_eval,
                                                int n,
                                                const double* x_nodes,
                                                const double* xi_nodes,
                                                const double* coeffs);


double calculate_parabolic_spline_discrepancy(double x_eval,
                                              double (*func)(double),
                                              int n,
                                              const double* x_nodes,
                                              const double* xi_nodes,
                                              const double* coeffs);

#endif // PARABOLIC_APPROXIMATION_H 
