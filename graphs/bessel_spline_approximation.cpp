#include <vector>
#include <cmath>
#include <limits>
#include <algorithm> // for std::upper_bound, std::min, std::max
#include <iostream>  // for std::cerr

#include "bessel_spline_approximation.h"

// Вспомогательная функция для вычисления производной d_i в узле x_i
// по трем точкам (x_prev, f_prev), (x_curr, f_curr), (x_next, f_next)
// Это формула Бесселя (19.4.7) стр. 102, переписанная для удобства
static double calculate_bessel_derivative(double x_prev, double f_prev,
                                          double x_curr, double f_curr,
                                          double x_next, double f_next) {
    double h_curr = x_curr - x_prev; // x_i - x_{i-1}
    double h_next = x_next - x_curr; // x_{i+1} - x_i

    if (std::abs(h_curr) < std::numeric_limits<double>::epsilon() ||
        std::abs(h_next) < std::numeric_limits<double>::epsilon() ||
        std::abs(x_next - x_prev) < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double term1_val = (f_next - f_curr) / h_next;
    double term2_val = (f_curr - f_prev) / h_curr;
    
    double derivative = (term1_val * h_curr + term2_val * h_next) / (h_curr + h_next);
    return derivative;
}


int make_bessel_spline_coefficients(int n, const double* x_nodes, const double* f_values, double* bessel_coeffs) {
    if (n < 1 || !x_nodes || !f_values || !bessel_coeffs) {
        // bessel_coeffs выделяется как 4*n в window.cpp, так что bessel_coeffs[0] безопасен если n >= 1
        if (bessel_coeffs && n >= 1) { 
             bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        }
        return -1;
    }

    if (n == 1) {
        // Для одной точки сплайн не строится в виде полиномов.
        // bessel_coeffs[0] (индекс 0 из массива размером 4*1=4) помечается NaN.
        bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN(); 
        return 0; 
    }

    // Проверка на монотонность узлов для n > 1 (т.е. n >= 2)
    for (int i = 0; i < n - 1; ++i) {
        if (x_nodes[i+1] - x_nodes[i] <= std::numeric_limits<double>::epsilon()) {
            // bessel_coeffs[0] безопасен, т.к. n >= 2
            bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    std::vector<double> d(n); // Массив для хранения производных d_i

    // 1. Вычисление производных d_i
    if (n == 2) { // Частный случай: 2 узла, 1 интервал. d_0 и d_1 - оба граничные.
        double x0_fict, f0_fict; 
        double x3_fict, f3_fict; 

        double h1 = x_nodes[1] - x_nodes[0];
        if (std::abs(h1) < std::numeric_limits<double>::epsilon()) { 
            // bessel_coeffs[0] безопасен, т.к. n = 2
            bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN(); 
            return -1;
        }

        x0_fict = x_nodes[0] - h1;
        f0_fict = f_values[0] - (f_values[1] - f_values[0]); 

        x3_fict = x_nodes[1] + h1;
        f3_fict = f_values[1] + (f_values[1] - f_values[0]); 

        d[0] = calculate_bessel_derivative(x0_fict, f0_fict, x_nodes[0], f_values[0], x_nodes[1], f_values[1]);
        d[1] = calculate_bessel_derivative(x_nodes[0], f_values[0], x_nodes[1], f_values[1], x3_fict, f3_fict);
        
        if(std::isnan(d[0]) || std::isnan(d[1])) {
             // bessel_coeffs[0] безопасен, т.к. n = 2
             bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN(); 
             return -1;
        }

    } else { // n >= 3
        double x_fict_left = x_nodes[0] - (x_nodes[1] - x_nodes[0]);
        double f_fict_left = f_values[0] - (f_values[1] - f_values[0]);
        d[0] = calculate_bessel_derivative(x_fict_left, f_fict_left, x_nodes[0], f_values[0], x_nodes[1], f_values[1]);

        double x_fict_right = x_nodes[n-1] + (x_nodes[n-1] - x_nodes[n-2]);
        double f_fict_right = f_values[n-1] + (f_values[n-1] - f_values[n-2]);
        d[n-1] = calculate_bessel_derivative(x_nodes[n-2], f_values[n-2], x_nodes[n-1], f_values[n-1], x_fict_right, f_fict_right);

        if(std::isnan(d[0]) || std::isnan(d[n-1])) {
            // bessel_coeffs[0] безопасен, т.к. n >= 3
            bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN(); 
            return -1;
        }

        for (int i = 1; i < n - 1; ++i) {
            d[i] = calculate_bessel_derivative(x_nodes[i-1], f_values[i-1], x_nodes[i], f_values[i], x_nodes[i+1], f_values[i+1]);
            if(std::isnan(d[i])) {
                // bessel_coeffs[0] безопасен, т.к. n >= 3
                bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN(); 
                return -1;
            }
        }
    }

    // 2. Вычисление коэффициентов кубических многочленов для каждого интервала [x_k, x_{k+1}]
    for (int k = 0; k < n - 1; ++k) { // k от 0 до n-2. Цикл выполняется если n >= 2.
        double xk = x_nodes[k];
        double xk1 = x_nodes[k+1];
        double fk = f_values[k];
        double fk1 = f_values[k+1];
        double dk_val = d[k]; // Используем dk_val, чтобы не конфликтовать с массивом d в С++20 (хотя здесь это не проблема)
        double dk1_val = d[k+1];

        double hk = xk1 - xk;
        if (std::abs(hk) < std::numeric_limits<double>::epsilon()) {
            // bessel_coeffs[0] безопасен, т.к. n >= 2
            bessel_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1; 
        }

        bessel_coeffs[4 * k + 0] = fk;                                  
        bessel_coeffs[4 * k + 1] = dk_val;                                  
        bessel_coeffs[4 * k + 2] = (3*(fk1-fk)/(hk) - 2*dk_val - dk1_val) / hk; 
        bessel_coeffs[4 * k + 3] = (dk_val + dk1_val - 2*(fk1-fk)/(hk)) / (hk*hk); 
    }
    
    return 0;
}


double calculate_bessel_spline_approximation(double x_eval, int n, const double* x_nodes, const double* bessel_coeffs) {
    if (n < 1 || !x_nodes || !bessel_coeffs) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (n == 1) {
        // Для n=1, make_bessel_spline_coefficients устанавливает bessel_coeffs[0] = NaN.
        // Аппроксимация одной точкой не является "сплайном" в обычном смысле.
        // Возвращаем NaN, чтобы указать, что метод Бесселя как таковой не построил полином.
        // Отрисовка в Window::paintEvent должна корректно обработать NaN (не рисовать).
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Проверка на флаг ошибки в коэффициентах (устанавливается в bessel_coeffs[0])
    // Этот доступ безопасен, т.к. bessel_coeffs выделяется как 4*n, и n >= 2 здесь.
    if (std::isnan(bessel_coeffs[0])) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto it = std::upper_bound(x_nodes, x_nodes + n, x_eval);
    int idx = std::distance(x_nodes, it); 

    int k; 
    if (idx == 0) { 
        k = 0;
    } else if (idx == n) { 
        k = n - 2;
    } else { 
        k = idx - 1;
    }
    k = std::max(0, std::min(k, n - 2));


    double xk_node = x_nodes[k]; // Используем xk_node, чтобы не конфликтовать с k из цикла
    double t = x_eval - xk_node;

    const double* c_ptr = &bessel_coeffs[4 * k]; 
    double approx_val = c_ptr[0] + t * (c_ptr[1] + t * (c_ptr[2] + t * c_ptr[3]));

    return approx_val;
}

double calculate_bessel_spline_discrepancy(double x_eval, func_t func, int n, const double* x_nodes, const double* bessel_coeffs) {
    if (!func) return std::numeric_limits<double>::quiet_NaN();
    double approx_val = calculate_bessel_spline_approximation(x_eval, n, x_nodes, bessel_coeffs);
    if (std::isnan(approx_val)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double true_val = func(x_eval);
    return std::abs(true_val - approx_val);
}
