#ifndef CHEBYSHEV_APPROXIMATION_H
#define CHEBYSHEV_APPROXIMATION_H

#include "functions.h" // Для func_t

// n - количество узлов (и коэффициентов alpha_0 ... alpha_{n-1})
// a, b - границы интервала аппроксимации
// original_function - указатель на исходную функцию f(x)
// chebyshev_coeffs - выходной массив для коэффициентов alpha_i (размер n)
// Возвращает 0 при успехе, -1 при ошибке (например, n < 1)
int make_chebyshev_coefficients(int n, double a, double b, func_t original_function, double* chebyshev_coeffs);

// x_eval - точка, в которой вычисляется значение
// n - количество коэффициентов
// a, b - границы интервала
// chebyshev_coeffs - массив коэффициентов (размер n)
// Возвращает значение аппроксимации или NaN при ошибке
double calculate_chebyshev_approximation(double x_eval, int n, double a, double b, const double* chebyshev_coeffs);

// x_eval - точка вычисления невязки
// func - исходная функция
// n, a, b, chebyshev_coeffs - параметры аппроксимации
// Возвращает абсолютную невязку |func(x_eval) - Approx(x_eval)| или NaN
double calculate_chebyshev_discrepancy(double x_eval, func_t func, int n, double a, double b, const double* chebyshev_coeffs);

#endif // CHEBYSHEV_APPROXIMATION_H 
