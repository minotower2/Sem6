#ifndef BESSEL_SPLINE_APPROXIMATION_H
#define BESSEL_SPLINE_APPROXIMATION_H

#include "functions.h" // Для func_t

// n - количество узлов (x_nodes, f_values)
// x_nodes - массив узлов x_0, ..., x_{n-1} (размер n, отсортирован)
// f_values - массив значений функции f(x_i) в узлах (размер n)
// bessel_coeffs - выходной массив для коэффициентов.
//                 Размер должен быть 4*(n-1), если n > 1.
//                 Для n=1 коэффициенты не вычисляются, bessel_coeffs[0] = NaN.
//                 Для каждого интервала [x_i, x_{i+1}], i=0..n-2:
//                 c[4*i+0]=a_i, c[4*i+1]=b_i, c[4*i+2]=c_i, c[4*i+3]=d_i
// Возвращает 0 при успехе, -1 при ошибке.
int make_bessel_spline_coefficients(int n, const double* x_nodes, const double* f_values, double* bessel_coeffs);

// x_eval - точка для вычисления значения
// n - количество исходных узлов
// x_nodes - массив узлов (размер n)
// bessel_coeffs - массив коэффициентов (размер 4*(n-1) для n > 1)
// Возвращает значение аппроксимации или NaN.
double calculate_bessel_spline_approximation(double x_eval, int n, const double* x_nodes, const double* bessel_coeffs);

// x_eval - точка для вычисления невязки
// func - исходная функция
// n, x_nodes, bessel_coeffs - параметры аппроксимации
// Возвращает абсолютную невязку или NaN.
double calculate_bessel_spline_discrepancy(double x_eval, func_t func, int n, const double* x_nodes, const double* bessel_coeffs);

#endif // BESSEL_SPLINE_APPROXIMATION_H 
