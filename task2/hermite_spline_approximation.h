#ifndef HERMITE_SPLINE_APPROXIMATION_H
#define HERMITE_SPLINE_APPROXIMATION_H


#include "functions.h"
// Построение коэффициентов кубического эрмитова сплайна по заданным узлам,
// значениям функции и значениям производных в узлах.
// Параметры:
//   n             - число узлов (n >= 1)
//   x_nodes       - массив узлов (строго возрастающие)
//   f_values      - массив значений функции в узлах
//   d_values      - массив значений производных в узлах
//   hermite_coeffs - выходной массив размером 4*(n-1) (для n=1 не используется)
// Возвращает 0 при успехе, -1 при ошибке (при ошибке hermite_coeffs[0] = NaN).
int make_cubic_hermite_coefficients(int n,
                                    const double* x_nodes,
                                    const double* f_values,
                                    const double* d_values,
                                    double* hermite_coeffs);

// Вычисление значения сплайна в произвольной точке.
double calculate_cubic_hermite_approximation(double x_eval,
                                             int n,
                                             const double* x_nodes,
                                             const double* hermite_coeffs);

// Вычисление абсолютной погрешности (разность между сплайном и заданной функцией func).
double calculate_cubic_hermite_discrepancy(double x_eval,
                                           func_t func,
                                           int n,
                                           const double* x_nodes,
                                           const double* hermite_coeffs);

#endif // HERMITE_SPLINE_APPROXIMATION_H
