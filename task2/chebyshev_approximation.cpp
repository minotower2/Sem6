#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate (хотя не используется напрямую)
#include <limits>  // for std::numeric_limits
#include <iostream> // for std::cerr (отладка)

#include "chebyshev_approximation.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Вспомогательная функция для вычисления значения многочлена Чебышева T_k(z) на [-1, 1]
// k - степень, z - точка в [-1, 1]
static double chebyshev_T(int k, double z) {
    if (k == 0) return 1.0;
    if (k == 1) return z;
    // Можно использовать std::cos(k * std::acos(z)), но для стабильности лучше рекуррентная формула
    // Однако, для прямого вычисления T_k(z_j) при построении g_ij, cos(k*acos(z)) может быть удобнее
    // так как z_j = cos(theta_j), то T_k(z_j) = cos(k*theta_j)
    if (z > 1.0) z = 1.0; // Защита для acos
    if (z < -1.0) z = -1.0;
    return std::cos(static_cast<double>(k) * std::acos(z));
}

// Вычисление коэффициентов разложения по многочленам Чебышева
// n_nodes - количество узлов (и коэффициентов alpha_0 ... alpha_{n_nodes-1}), соответствует N в книге
// a, b - границы интервала
// original_function - указатель на исходную функцию f(x)
// chebyshev_coeffs - выходной массив для коэффициентов alpha_i (размер n_nodes)
int make_chebyshev_coefficients(int n_nodes, double a, double b, func_t original_function, double* chebyshev_coeffs) {
    if (n_nodes < 1) {
        if (chebyshev_coeffs && n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    if (!chebyshev_coeffs || !original_function) {
        // Если chebyshev_coeffs != nullptr, но n_nodes < 1, то первый элемент уже помечен NaN.
        // Если chebyshev_coeffs == nullptr, но n_nodes >= 1, это ошибка.
        // Пометим, если возможно.
        if (chebyshev_coeffs && n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    if (std::abs(b - a) < std::numeric_limits<double>::epsilon()) { // Отрезок нулевой длины
         if (n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }


    std::vector<double> x_cheb_nodes(n_nodes); // Чебышевские узлы на [a,b]
    std::vector<double> f_at_nodes(n_nodes);   // Значения f(x) в этих узлах
    std::vector<double> z_cheb_args(n_nodes);  // Аргументы для T_k на [-1,1]

    // Вычисляем чебышевские узлы x_j и значения z_j = cos( (2j-1)pi / (2N) )
    // Индексация узлов в книге x_m, m=1..N. В C++ j=0..N-1.
    // x_m = (a+b)/2 + (b-a)/2 * cos( (2m-1)pi / (2N) )
    // z_m = cos( (2m-1)pi / (2N) )
    for (int j = 0; j < n_nodes; ++j) {
        double cos_arg = (2.0 * (j + 1.0) - 1.0) * M_PI / (2.0 * n_nodes);
        z_cheb_args[j] = std::cos(cos_arg);
        x_cheb_nodes[j] = 0.5 * (a + b) + 0.5 * (b - a) * z_cheb_args[j];
        f_at_nodes[j] = original_function(x_cheb_nodes[j]);
    }

    // Вычисляем коэффициенты alpha_i (в коде chebyshev_coeffs[i])
    // alpha_i = (norm_factor / n_nodes) * sum_{j=0}^{n_nodes-1} f_at_nodes[j] * T_i(z_cheb_args[j])
    // norm_factor = 1 для i=0, 2 для i > 0

    for (int i = 0; i < n_nodes; ++i) { // Индекс коэффициента alpha_i
        double current_sum = 0;
        for (int j = 0; j < n_nodes; ++j) { // Суммирование по узлам x_j (в C++ индексы 0..n_nodes-1)
            // T_i(z_cheb_args[j]) = T_i(cos( (2(j+1)-1)pi / (2N) ))
            // T_i(cos(theta)) = cos(i*theta)
            double cos_theta_j = z_cheb_args[j]; // Это cos( (2(j+1)-1)pi / (2N) )
            double T_i_val_at_z_j = chebyshev_T(i, cos_theta_j); // T_i(z_j)
            current_sum += f_at_nodes[j] * T_i_val_at_z_j;
        }
        if (i == 0) {
            chebyshev_coeffs[i] = current_sum / static_cast<double>(n_nodes);
        } else {
            chebyshev_coeffs[i] = 2.0 * current_sum / static_cast<double>(n_nodes);
        }
    }
    return 0;
}

// Вычисление значения аппроксимации по многочленам Чебышева
// x_eval - точка, в которой вычисляется значение
// n_coeffs - количество коэффициентов (степень многочлена N-1 = n_coeffs-1)
// a, b - границы интервала
// chebyshev_coeffs - массив коэффициентов alpha_i (размер n_coeffs)
double calculate_chebyshev_approximation(double x_eval, int n_coeffs, double a, double b, const double* chebyshev_coeffs) {
    if (n_coeffs < 1 || !chebyshev_coeffs) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(b - a) < std::numeric_limits<double>::epsilon()) { // Отрезок нулевой длины
        if (n_coeffs == 1) return chebyshev_coeffs[0]; // Константа
        return std::numeric_limits<double>::quiet_NaN(); // Неопределенность
    }
    if (std::isnan(chebyshev_coeffs[0])) { // Проверка на флаг ошибки от make_
        return std::numeric_limits<double>::quiet_NaN();
    }


    double z_eval = (2.0 * x_eval - (a + b)) / (b - a);
    // Ограничиваем z_eval диапазоном [-1, 1] для стабильности T_k(z) при экстраполяции
    // Хотя для полиномов это не строго обязательно, но T_k(z) через acos этого требует
    if (z_eval > 1.0 && z_eval < 1.0 + 1e-9) z_eval = 1.0; // Небольшая коррекция для чисел на границе
    if (z_eval < -1.0 && z_eval > -1.0 - 1e-9) z_eval = -1.0;


    // Схема Кленшоу для S = sum_{k=0}^{N-1} alpha_k T_k(z)
    // N = n_coeffs
    if (n_coeffs == 1) {
        return chebyshev_coeffs[0]; // T_0(z) = 1
    }

    double y_k_plus_2 = 0.0;
    double y_k_plus_1 = 0.0; // y_N = 0
    
    // Если N-1 это индекс последнего коэффициента, то N-1 = n_coeffs - 1
    // y_{n_coeffs} = 0
    // y_{n_coeffs-1} = chebyshev_coeffs[n_coeffs-1]
    
    y_k_plus_1 = chebyshev_coeffs[n_coeffs - 1]; // y_{N-1}

    if (n_coeffs == 2) { // alpha_0 T_0 + alpha_1 T_1 = alpha_0 + alpha_1 * z
        // y_2 = 0
        // y_1 = alpha_1
        // y_0 = alpha_0 + 2*z*y_1 - y_2 = alpha_0 + 2*z*alpha_1
        // Result = y_0 - z*y_1 = alpha_0 + 2*z*alpha_1 - z*alpha_1 = alpha_0 + z*alpha_1
        return chebyshev_coeffs[0] + chebyshev_coeffs[1] * z_eval;
    }


    // Цикл для k = N-2, ..., 0
    // В C++ k_idx = n_coeffs-2 ... 0
    for (int k_idx = n_coeffs - 2; k_idx >= 0; --k_idx) {
        double y_k = chebyshev_coeffs[k_idx] + 2.0 * z_eval * y_k_plus_1 - y_k_plus_2;
        y_k_plus_2 = y_k_plus_1;
        y_k_plus_1 = y_k;
    }
    // После цикла: y_k_plus_1 = y_0, y_k_plus_2 = y_1
    // Result = y_0 - z*y_1
    return y_k_plus_1 - z_eval * y_k_plus_2;
}

// Вычисление невязки для аппроксимации Чебышева
double calculate_chebyshev_discrepancy(double x_eval, func_t func, int n, double a, double b, const double* chebyshev_coeffs) {
    if (!func) return std::numeric_limits<double>::quiet_NaN();
    double approx_val = calculate_chebyshev_approximation(x_eval, n, a, b, chebyshev_coeffs);
    if (std::isnan(approx_val)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double true_val = func(x_eval);
    return std::abs(true_val - approx_val);
} 
