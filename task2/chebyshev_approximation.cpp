#include <vector>
#include <cmath>
#include <numeric> 
#include <limits>  
#include <iostream> 

#include "chebyshev_approximation.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double chebyshev_T(int k, double z) {
    if (k == 0) return 1.0;
    if (k == 1) return z;
    if (z > 1.0) z = 1.0; // Защита для acos
    if (z < -1.0) z = -1.0;
    return std::cos(static_cast<double>(k) * std::acos(z));
}

int make_chebyshev_coefficients(int n_nodes, double a, double b, func_t original_function, double* chebyshev_coeffs) {
    if (n_nodes < 1) {
        if (chebyshev_coeffs && n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    if (!chebyshev_coeffs || !original_function) {
        if (chebyshev_coeffs && n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    if (std::abs(b - a) < std::numeric_limits<double>::epsilon()) { // Отрезок нулевой длины
         if (n_nodes > 0) chebyshev_coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }


    std::vector<double> x_cheb_nodes(n_nodes); 
    std::vector<double> f_at_nodes(n_nodes);   
    std::vector<double> z_cheb_args(n_nodes);  

    for (int j = 0; j < n_nodes; ++j) {
        double cos_arg = (2.0 * (j + 1.0) - 1.0) * M_PI / (2.0 * n_nodes);
        z_cheb_args[j] = std::cos(cos_arg);
        x_cheb_nodes[j] = 0.5 * (a + b) + 0.5 * (b - a) * z_cheb_args[j];
        f_at_nodes[j] = original_function(x_cheb_nodes[j]);
    }


    for (int i = 0; i < n_nodes; ++i) { 
        double current_sum = 0;
        for (int j = 0; j < n_nodes; ++j) { 
            double cos_theta_j = z_cheb_args[j]; 
            double T_i_val_at_z_j = chebyshev_T(i, cos_theta_j); 
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

double calculate_chebyshev_approximation(double x_eval, int n_coeffs, double a, double b, const double* chebyshev_coeffs) {
    if (n_coeffs < 1 || !chebyshev_coeffs) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(b - a) < std::numeric_limits<double>::epsilon()) { 
        if (n_coeffs == 1) return chebyshev_coeffs[0]; 
        return std::numeric_limits<double>::quiet_NaN(); 
    }
    if (std::isnan(chebyshev_coeffs[0])) { 
        return std::numeric_limits<double>::quiet_NaN();
    }


    double z_eval = (2.0 * x_eval - (a + b)) / (b - a);
    if (z_eval > 1.0 && z_eval < 1.0 + 1e-9) z_eval = 1.0; 
    if (z_eval < -1.0 && z_eval > -1.0 - 1e-9) z_eval = -1.0;


    if (n_coeffs == 1) {
        return chebyshev_coeffs[0]; 
    }

    double y_k_plus_2 = 0.0;
    double y_k_plus_1 = 0.0; 
    
    
    y_k_plus_1 = chebyshev_coeffs[n_coeffs - 1]; // y_{N-1}

    if (n_coeffs == 2) { 
        return chebyshev_coeffs[0] + chebyshev_coeffs[1] * z_eval;
    }


    for (int k_idx = n_coeffs - 2; k_idx >= 0; --k_idx) {
        double y_k = chebyshev_coeffs[k_idx] + 2.0 * z_eval * y_k_plus_1 - y_k_plus_2;
        y_k_plus_2 = y_k_plus_1;
        y_k_plus_1 = y_k;
    }
    return y_k_plus_1 - z_eval * y_k_plus_2;
}

double calculate_chebyshev_discrepancy(double x_eval, func_t func, int n, double a, double b, const double* chebyshev_coeffs) {
    if (!func) return std::numeric_limits<double>::quiet_NaN();
    double approx_val = calculate_chebyshev_approximation(x_eval, n, a, b, chebyshev_coeffs);
    if (std::isnan(approx_val)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double true_val = func(x_eval);
    return std::abs(true_val - approx_val);
} 
