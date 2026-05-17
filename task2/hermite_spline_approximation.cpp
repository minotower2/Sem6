#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

#include "hermite_spline_approximation.h"

int make_cubic_hermite_coefficients(int n, 
                                    const double* x_nodes, 
                                    const double* f_values, 
                                    const double* d_values,   // производные в узлах
                                    double* hermite_coeffs) {
    if (n < 1 || !x_nodes || !f_values || !d_values || !hermite_coeffs) {
        if (hermite_coeffs && n >= 1)
            hermite_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    if (n == 1) {
        // Для одного узла сплайн не определён (нужен хотя бы один интервал)
        hermite_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return 0; // не ошибка, но аппроксимация невозможна
    }

    // Проверка строгого возрастания узлов
    for (int i = 0; i < n - 1; ++i) {
        if (x_nodes[i+1] - x_nodes[i] <= std::numeric_limits<double>::epsilon()) {
            hermite_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Для каждого интервала [x_k, x_{k+1}] строим кубический полином Эрмита
    for (int k = 0; k < n - 1; ++k) {
        double xk  = x_nodes[k];
        double xk1 = x_nodes[k+1];
        double fk  = f_values[k];
        double fk1 = f_values[k+1];
        double dk  = d_values[k];
        double dk1 = d_values[k+1];
        double h   = xk1 - xk;

        if (std::abs(h) < std::numeric_limits<double>::epsilon()) {
            hermite_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }

        // Коэффициенты кубического полинома: p(t) = c0 + c1*t + c2*t^2 + c3*t^3, t = x - xk
        hermite_coeffs[4*k + 0] = fk;                                 // c0
        hermite_coeffs[4*k + 1] = dk;                                 // c1
        hermite_coeffs[4*k + 2] = (3.0*(fk1 - fk)/h - 2.0*dk - dk1) / h;
        hermite_coeffs[4*k + 3] = (dk + dk1 - 2.0*(fk1 - fk)/h) / (h*h);
    }

    return 0;
}

double calculate_cubic_hermite_approximation(double x_eval, 
                                             int n, 
                                             const double* x_nodes, 
                                             const double* hermite_coeffs) {
    if (n < 1 || !x_nodes || !hermite_coeffs)
        return std::numeric_limits<double>::quiet_NaN();

    if (n == 1)
        return std::numeric_limits<double>::quiet_NaN();

    // Проверка на ошибку в коэффициентах (hermite_coeffs[0] == NaN)
    if (std::isnan(hermite_coeffs[0]))
        return std::numeric_limits<double>::quiet_NaN();

    // Поиск интервала, содержащего x_eval
    auto it = std::upper_bound(x_nodes, x_nodes + n, x_eval);
    int idx = std::distance(x_nodes, it);
    int k;
    if (idx == 0) {
        k = 0;                       // экстраполяция влево
    } else if (idx == n) {
        k = n - 2;                   // экстраполяция вправо
    } else {
        k = idx - 1;                 // обычный интервал
    }
    k = std::max(0, std::min(k, n - 2));

    double t = x_eval - x_nodes[k];
    const double* c = &hermite_coeffs[4*k];
    return c[0] + t * (c[1] + t * (c[2] + t * c[3]));
}

double calculate_cubic_hermite_discrepancy(double x_eval, 
                                           func_t func, 
                                           int n, 
                                           const double* x_nodes, 
                                           const double* hermite_coeffs) {
    if (!func)
        return std::numeric_limits<double>::quiet_NaN();

    double approx = calculate_cubic_hermite_approximation(x_eval, n, x_nodes, hermite_coeffs);
    if (std::isnan(approx))
        return std::numeric_limits<double>::quiet_NaN();

    double exact = func(x_eval);
    return std::abs(exact - approx);
}
