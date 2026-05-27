#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

static bool solve_tridiagonal(int n, const double* a, const double* b, const double* c,
                              const double* d, double* x) {
    // a, b, c, d – массивы длины n, a[0] и c[n-1] не используются
    std::vector<double> alpha(n), beta(n);
    alpha[0] = -c[0] / b[0];
    beta[0]  =  d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * alpha[i-1];
        if (std::abs(denom) < std::numeric_limits<double>::epsilon())
            return false;
        alpha[i] = -c[i] / denom;
        beta[i]  = (d[i] - a[i] * beta[i-1]) / denom;
    }
    x[n-1] = beta[n-1];
    for (int i = n-2; i >= 0; --i)
        x[i] = alpha[i] * x[i+1] + beta[i];
    return true;
}

int make_parabolic_spline_coefficients(int n,
                                       const double* x_nodes,
                                       const double* f_values,
                                       double* xi_nodes,
                                       double* coeffs) {
    if (n < 2 || !x_nodes || !f_values || !xi_nodes || !coeffs) {
        if (coeffs && n >= 1) coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Проверка строгого возрастания узлов
    for (int i = 0; i < n-1; ++i) {
        if (x_nodes[i+1] - x_nodes[i] <= std::numeric_limits<double>::epsilon()) {
            coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // 1. Построение дополнительных узлов ξ
    double h_left  = x_nodes[1] - x_nodes[0];
    double h_right = x_nodes[n-1] - x_nodes[n-2];
    xi_nodes[0] = x_nodes[0] - 0.5 * h_left;          // ξ₁
    xi_nodes[n] = x_nodes[n-1] + 0.5 * h_right;       // ξ_{n+1}
    for (int i = 1; i < n; ++i) {
        xi_nodes[i] = 0.5 * (x_nodes[i-1] + x_nodes[i]);  // ξ_{i+1}
    }
    // После этого xi_nodes[0]=ξ₁, xi_nodes[1]=ξ₂, ..., xi_nodes[n-1]=ξ_n, xi_nodes[n]=ξ_{n+1}

    // 2. Подготовка системы для v_i (размер N = n+1)
    int N = n + 1;
    std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), d(N, 0.0);

    // Вспомогательные массивы для длин
    std::vector<double> h_i(n), z_i(n); // h_i = x_i - ξ_i, z_i = ξ_{i+1} - x_i, i=1..n
    for (int i = 0; i < n; ++i) {
        h_i[i] = x_nodes[i] - xi_nodes[i];
        z_i[i] = xi_nodes[i+1] - xi_nodes[i];
        if (std::abs(h_i[i]) < std::numeric_limits<double>::epsilon() ||
            std::abs(z_i[i]) < std::numeric_limits<double>::epsilon()) {
            coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Граничное условие при i=1 (P₁''(x₁)=0)
    // Уравнение: v₁/(x₁-ξ₁) + v₂/(ξ₂-x₁) = f(x₁)*(1/(x₁-ξ₁) + 1/(ξ₂-x₁))
    double denom1 = 1.0/h_i[0] + 1.0/z_i[0];
    b[0] = 1.0/h_i[0];      // коэффициент при v₁
    c[0] = 1.0/z_i[0];      // коэффициент при v₂
    d[0] = f_values[0] * denom1;
    a[0] = 0.0;

    // Внутренние уравнения для i=2..n (условия непрерывности производной в ξ_i)
    for (int i = 1; i < n; ++i) {   // i соответствует номеру ξ_{i+1} в книге, индекс i в C++ от 1 до n-1
        double h_left_span = x_nodes[i-i] - xi_nodes[i-1];   // x_{i} - ξ_{i}
        double z_left_span = xi_nodes[i] - x_nodes[i-1];    // ξ_{i+1} - x_{i}
        double h_right_span = x_nodes[i] - xi_nodes[i];  // x_{i+1} - ξ_{i+1}
        double z_right_span = xi_nodes[i+1] - x_nodes[i]; // ξ_{i+2} - x_{i+1}

        if (std::abs(h_left_span) < std::numeric_limits<double>::epsilon() ||
            std::abs(z_left_span) < std::numeric_limits<double>::epsilon() ||
            std::abs(h_right_span) < std::numeric_limits<double>::epsilon() ||
            std::abs(z_right_span) < std::numeric_limits<double>::epsilon()) {
            coeffs[0] = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }

        a[i] = 1.0/h_left_span - 1.0/(xi_nodes[i] - xi_nodes[i-1]);
        b[i]  = 1.0/h_right_span + 1.0/(xi_nodes[i+1] - xi_nodes[i]);
        c[i] = 1.0/z_right_span + 1.0/(xi_nodes[i+1] - xi_nodes[i]);
        d[i] = f_values[i-1] * (1.0/h_left_span + 1.0/z_left_span) +
                     f_values[i]   * (1.0/h_right_span + 1.0/z_right_span);

    }

    // Граничное условие при i=n (P_n''(x_n)=0)
    // Уравнение: v_n/(x_n-ξ_n) + v_{n+1}/(ξ_{n+1}-x_n) = f(x_n)*(1/(x_n-ξ_n) + 1/(ξ_{n+1}-x_n))
    double denom_last = 1.0/h_i[n-1] + 1.0/z_i[n-1];
    a[n] = 1.0/h_i[n-1];      // коэффициент при v_n
    b[n] = 1.0/z_i[n-1];      // коэффициент при v_{n+1}
    d[n] = f_values[n-1] * denom_last;
    c[n] = 0.0;

    // 3. Решение трёхдиагональной системы для v (размер N)
    std::vector<double> v(N);
    if (!solve_tridiagonal(N, a.data(), b.data(), c.data(), d.data(), v.data())) {
        coeffs[0] = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // 4. Вычисление коэффициентов квадратичных полиномов для каждого интервала [ξ_i, ξ_{i+1}]
    for (int i = 0; i < n; ++i) {
        double x = x_nodes[i];
        double h = x - xi_nodes[i];          // h_i
        double z = xi_nodes[i+1] - x;     // z_i

        double inv_h = 1.0 / h;
        double inv_z = 1.0 / z;
        double inv_hz = 1.0 / (h + z);

        // Коэффициенты квадратичного полинома: P(x) = c0 + c1*t + c2*t^2, t = x - ξ_i

        coeffs[3*i + 0] = v[i];
        coeffs[3*i + 1] = -v[i] * (inv_h + inv_hz) + f_values[i] * (inv_h + inv_z) - v[i+1] * (inv_z - inv_hz);
        coeffs[3*i + 2] = v[i] * inv_hz - f_values[i] * inv_hz * (inv_h + inv_z) + v[i+1] * inv_hz * inv_z;
    }

    return 0;
}

// Вычисление значения параболического сплайна в точке x_eval
double calculate_parabolic_spline_approximation(double x_eval,
                                                int n,
                                                const double* x_nodes,
                                                const double* xi_nodes,
                                                const double* coeffs) {
    if (n < 2 || !x_nodes || !xi_nodes || !coeffs)
        return std::numeric_limits<double>::quiet_NaN();

    // Поиск интервала [ξ_i, ξ_{i+1}], содержащего x_eval
    // Используем бинарный поиск по xi_nodes (длина n+1)
    auto it = std::upper_bound(xi_nodes, xi_nodes + n + 1, x_eval);
    int idx = std::distance(xi_nodes, it);
    int i;
    if (idx == 0) {
        i = 0;                     // экстраполяция влево
    } else if (idx == n+1) {
        i = n-1;                   // экстраполяция вправо
    } else {
        i = idx - 1;
    }
    i = std::max(0, std::min(i, n-1));

    double t = x_eval - xi_nodes[i];
    const double* c = &coeffs[3*i];
    return c[0] + t * (c[1] + t * c[2]);
}

// Вычисление невязки (абсолютной погрешности) между сплайном и функцией func в точке
double calculate_parabolic_spline_discrepancy(double x_eval,
                                              double (*func)(double),
                                              int n,
                                              const double* x_nodes,
                                              const double* xi_nodes,
                                              const double* coeffs) {
    if (!func) return std::numeric_limits<double>::quiet_NaN();
    double approx = calculate_parabolic_spline_approximation(x_eval, n, x_nodes, xi_nodes, coeffs);
    if (std::isnan(approx)) return std::numeric_limits<double>::quiet_NaN();
    double exact = func(x_eval);
    return std::abs(exact - approx);
}
