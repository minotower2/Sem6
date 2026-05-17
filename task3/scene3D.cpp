#include "scene3D.h"
#include <QtGui>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring> // For std::memset if used, or memcpy
#include <algorithm>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UNUSED(x) (void) x
#define EPS_SMALL 1e-9     // Малая величина для избежания деления на ноль и сравнения с нулем

// --- Тестовые функции (только сами функции, производные не нужны для этого метода сплайнов) ---
double f0_func(double x, double y) { UNUSED(x); UNUSED(y); return 1.0; }
double f1_func(double x, double y) { UNUSED(y); return x; }
double f2_func(double x, double y) { UNUSED(x); return y; }
double f3_func(double x, double y) { return x + y; }
double f4_func(double x, double y) { return sqrt(x*x + y*y); }
double f5_func(double x, double y) { return x*x + y*y; }
double f6_func(double x, double y) { return exp(x*x - y*y); }
double f7_func(double x, double y) { return 1.0 / (25.0*(x*x + y*y) + 1.0); }


Scene3D::Scene3D(char *argv[], QWidget *parent)
    : QGLWidget(parent), xRot(-90), yRot(0), zRot(0), zTra(0), nSca(1),
      spline_x_nodes(nullptr), spline_y_nodes(nullptr),
      f_val(nullptr), fx_val(nullptr), fy_val(nullptr), fxy_val(nullptr)
{
    nx_param = 10; ny_param = 10; id = 0; eps_param = 1e-5;
    a = -0.5; b = 0.5; c = -0.5; d = 0.5;

    int provided_argc = 0;
    for(int i=0; argv[i] != nullptr; ++i) provided_argc++;

    if (parse_command_line(provided_argc, argv)) {
        std::cerr << "Warning: Command line parsing failed or wrong arguments, using defaults." << std::endl;
        // Ensure nx/ny_param are at least MIN_SPLINE_NODES if parsing fails with smaller values
        if (nx_param < MIN_SPLINE_NODES) nx_param = MIN_SPLINE_NODES;
        if (ny_param < MIN_SPLINE_NODES) ny_param = MIN_SPLINE_NODES;
    }
    
    current_nx = nx_param;
    current_ny = ny_param;

    p = 0; s = 0; max_metod = 0; repeat = 0; scale = 0;

    change_func(); // This will call computeSplineData
}

Scene3D::~Scene3D()
{
    cleanupSplineData();
}

void Scene3D::allocateSplineData(int rows_nx, int cols_ny) {
    cleanupSplineData();

    current_nx = rows_nx;
    current_ny = cols_ny;

    spline_x_nodes = new double[current_nx];
    spline_y_nodes = new double[current_ny];

    f_val = new double*[current_nx];
    fx_val = new double*[current_nx];
    fy_val = new double*[current_nx];
    fxy_val = new double*[current_nx];

    for (int i = 0; i < current_nx; ++i) {
        f_val[i] = new double[current_ny]();
        fx_val[i] = new double[current_ny]();
        fy_val[i] = new double[current_ny]();
        fxy_val[i] = new double[current_ny]();
    }
}

void Scene3D::cleanupSplineData()
{
    if (f_val) {
        for (int i = 0; i < current_nx; ++i) delete[] f_val[i];
        delete[] f_val; f_val = nullptr;
    }
    if (fx_val) {
        for (int i = 0; i < current_nx; ++i) delete[] fx_val[i];
        delete[] fx_val; fx_val = nullptr;
    }
    if (fy_val) {
        for (int i = 0; i < current_nx; ++i) delete[] fy_val[i];
        delete[] fy_val; fy_val = nullptr;
    }
    if (fxy_val) {
        for (int i = 0; i < current_nx; ++i) delete[] fxy_val[i];
        delete[] fxy_val; fxy_val = nullptr;
    }
    delete[] spline_x_nodes; spline_x_nodes = nullptr;
    delete[] spline_y_nodes; spline_y_nodes = nullptr;
    current_nx = 0;
    current_ny = 0;
}

void Scene3D::change_func()
{
    id = (id + 1) % 8;
    switch (id) {
        case 0: f_name = "f(x,y) = 1"; f = f0_func; break;
        case 1: f_name = "f(x,y) = x"; f = f1_func; break;
        case 2: f_name = "f(x,y) = y"; f = f2_func; break;
        case 3: f_name = "f(x,y) = x+y"; f = f3_func; break;
        case 4: f_name = "f(x,y) = sqrt(x^2+y^2)"; f = f4_func; break;
        case 5: f_name = "f(x,y) = x^2+y^2"; f = f5_func; break;
        case 6: f_name = "f(x,y) = exp(x^2-y^2)"; f = f6_func; break;
        case 7: f_name = "f(x,y) = 1/(25*(x^2+y^2)+1)"; f = f7_func; break;
        default:f = f0_func; break;
    }
    fmax_orig = max_original_func();
    computeSplineData();
    updateGL();
}

void Scene3D::change_graph() { gr = (gr + 1) % 3; updateGL(); }

void Scene3D::increase_param_scale() {
    a /= 2.0; b /= 2.0; c /= 2.0; d /= 2.0; scale--;
    fmax_orig = max_original_func(); computeSplineData(); updateGL();
}
void Scene3D::decrease_param_scale() {
    a *= 2.0; b *= 2.0; c *= 2.0; d *= 2.0; scale++;
    fmax_orig = max_original_func(); computeSplineData(); updateGL();
}

void Scene3D::increase_nx_ny_param() {
    if ((nx_param <= 3000) && (ny_param <= 3000)){
        nx_param *= 2; ny_param *= 2; s--;
    }
    computeSplineData(); updateGL();
}
void Scene3D::decrease_nx_ny_param() {
    int new_nx = nx_param / 2;
    int new_ny = ny_param / 2;
    if (new_nx >= MIN_SPLINE_NODES && new_ny >= MIN_SPLINE_NODES) {
        nx_param = new_nx; ny_param = new_ny; s++;
    } else if (nx_param > MIN_SPLINE_NODES && ny_param > MIN_SPLINE_NODES) {
      // Allow decreasing if one dim is still large enough, but the other becomes MIN_SPLINE_NODES
      if (new_nx < MIN_SPLINE_NODES && nx_param > MIN_SPLINE_NODES) nx_param = MIN_SPLINE_NODES; else if (new_nx >= MIN_SPLINE_NODES) nx_param = new_nx;
      if (new_ny < MIN_SPLINE_NODES && ny_param > MIN_SPLINE_NODES) ny_param = MIN_SPLINE_NODES; else if (new_ny >= MIN_SPLINE_NODES) ny_param = new_ny;
      if (nx_param != current_nx || ny_param != current_ny) s++; // if changed
    }


    computeSplineData(); updateGL();
}

void Scene3D::increase_error_p() { p++; computeSplineData(); updateGL(); }
void Scene3D::decrease_error_p() { p--; computeSplineData(); updateGL(); }
void Scene3D::plus_angle(){ p_angle(); updateGL(); }
void Scene3D::minus_angle(){ m_angle(); updateGL(); }

int Scene3D::parse_command_line(int argc, char **argv) {
    if (argc < 6) return -1;
    FILE *file_ptr = fopen(argv[1], "r");
    if (!file_ptr) return -1;

    char line_buffer[100];
    double parsed_coords[4];
    int coords_count = 0;
    while (fgets(line_buffer, sizeof(line_buffer), file_ptr) != NULL) {
        if (line_buffer[0] == '#' || line_buffer[0] == ' ' || line_buffer[0] == '\n' || strlen(line_buffer) < 2) continue;
        std::istringstream iss(line_buffer);
        double tmp_coord1, tmp_coord2;
        if (iss >> tmp_coord1 >> tmp_coord2 && coords_count < 4) {
            parsed_coords[coords_count++] = tmp_coord1;
            parsed_coords[coords_count++] = tmp_coord2;
        }
        if (coords_count == 4) break;
    }
    fclose(file_ptr);
    if (coords_count != 4) return -1;

    this->a = parsed_coords[0]; this->d = parsed_coords[1];
    this->b = parsed_coords[2]; this->c = parsed_coords[3];
    if (this->a > this->b) std::swap(this->a, this->b);
    if (this->c > this->d) std::swap(this->c, this->d);

    char *eps_end_ptr;
    eps_param = strtod(argv[5], &eps_end_ptr);

    if (sscanf(argv[2], "%d", &nx_param) != 1 ||
        sscanf(argv[3], "%d", &ny_param) != 1 ||
        nx_param < MIN_SPLINE_NODES || ny_param < MIN_SPLINE_NODES  || // Check against MIN_SPLINE_NODES
        sscanf(argv[4], "%d", &id) != 1 ||
        eps_param <= 0 || id < 0 || id > 7 || eps_end_ptr == argv[5]) {
        return -2;
    }
    current_nx = nx_param;
    current_ny = ny_param;
    return 0;
}

double Scene3D::max_original_func() {
    if (!f) return 0.0;
    double val, max_val = 0.0;
    int eval_points_x = std::max(100, current_nx * 2);
    int eval_points_y = std::max(100, current_ny * 2);
    double dx = (b - a) / (eval_points_x -1 < 1 ? 1 : eval_points_x -1);
    double dy = (d - c) / (eval_points_y -1 < 1 ? 1 : eval_points_y -1);

    for (int i = 0; i < eval_points_x; ++i) {
        double cur_x = a + i * dx;
        for (int j = 0; j < eval_points_y; ++j) {
            double cur_y = c + j * dy;
            val = fabs(f(cur_x, cur_y));
            if (val > max_val) max_val = val;
        }
    }
    if (eval_points_x == 1 && eval_points_y == 1) max_val = fabs(f(a,c));
    return max_val;
}


// Решает СЛАУ Ax=B для трехдиагональной матрицы методом прогонки.
// a - поддиагональ (a[i] для i-го уравнения, i=1..N-1, a[0] не используется),
// d_mod - главная диагональ (d_mod[i] для i-го уравнения, i=0..N-1), МОДИФИЦИРУЕТСЯ
// c - наддиагональ (c[i] для i-го уравнения, i=0..N-2),
// b_mod - правая часть, МОДИФИЦИРУЕТСЯ
// x_sol - решение (размер N_sys)
// N_sys - размер системы.
bool Scene3D::solveTridiagonalSystem(const double* a, double* d_mod, const double* c, double* b_mod,
                                     double* x_sol, int N_sys)
{
    if (N_sys <= 0) return false;
    if (N_sys == 1) {
        if (std::abs(d_mod[0]) < EPS_SMALL) return false;
        x_sol[0] = b_mod[0] / d_mod[0];
        return true;
    }

    // Прямой ход (модифицируем c (как alpha) и b_mod (как beta))
    // c_prime (alpha) и b_prime (beta) хранятся на месте c и b_mod для экономии памяти
    // В стандартных обозначениях метода прогонки:
    // alpha[0] = c[0]/d[0], beta[0] = b[0]/d[0]
    // alpha[i] = c[i]/(d[i] - a[i-1]*alpha[i-1])
    // beta[i]  = (b[i] - a[i-1]*beta[i-1]) / (d[i] - a[i-1]*alpha[i-1])
    // Но здесь индексация массивов a,d,c немного другая, адаптируем:
    // a_i x_{i-1} + d_i x_i + c_i x_{i+1} = b_i (для уравнений i=0..N_sys-1)
    // a[0] не существует, c[N_sys-1] не существует.
    // Переменные для прогоночных коэффициентов (будут храниться в c и b_mod)
    double* c_prime = new double[N_sys]; // alpha
    double* b_prime = new double[N_sys]; // beta
    
    if (std::abs(d_mod[0]) < EPS_SMALL) { delete[] c_prime; delete[] b_prime; return false; }
    c_prime[0] = c[0] / d_mod[0];
    b_prime[0] = b_mod[0] / d_mod[0];

    for (int i = 1; i < N_sys; ++i) {
        double den = d_mod[i] - a[i-1] * c_prime[i-1]; // a[i-1] - это `a_i` из стандартной формулы прогонки для i-го уравнения
        if (std::abs(den) < EPS_SMALL) { delete[] c_prime; delete[] b_prime; return false; }
        if (i < N_sys - 1) { // c[i] существует
             c_prime[i] = c[i] / den;
        } else { // для последнего уравнения c_prime не нужен
             c_prime[i] = 0; // или не использовать
        }
        b_prime[i] = (b_mod[i] - a[i-1] * b_prime[i-1]) / den;
    }

    // Обратный ход
    x_sol[N_sys-1] = b_prime[N_sys-1];
    for (int i = N_sys - 2; i >= 0; --i) {
        x_sol[i] = b_prime[i] - c_prime[i] * x_sol[i+1];
    }
    delete[] c_prime;
    delete[] b_prime;
    return true;
}


// Экстраполирует вторую производную M на границе (M0 или M_{N-1}) используя 3 точки.
// coords: x_0, x_1, ..., x_{N_total-1}
// values: f_0, f_1, ..., f_{N_total-1}
// N_total: общее количество точек в 1D наборе
// at_start_node: true для M_0 (использует точки 0,1,2), false для M_{N_total-1} (использует N_total-3, N_total-2, N_total-1)
double Scene3D::extrapolate_M_boundary(const double* coords, const double* values, int N_total, bool at_start_node) {
    if (N_total < 3) { // Недостаточно точек для 3-точечной экстраполяции
        return 0.0;    // Условие "естественного" сплайна
    }

    double x0, x1, x2;
    double f0, f1, f2;

    if (at_start_node) { // Для M_0 используем точки x_0, x_1, x_2
        x0 = coords[0]; x1 = coords[1]; x2 = coords[2];
        f0 = values[0]; f1 = values[1]; f2 = values[2];
    } else { // Для M_{N-1} используем точки x_{N-3}, x_{N-2}, x_{N-1}
        x0 = coords[N_total - 3]; x1 = coords[N_total - 2]; x2 = coords[N_total - 1];
        f0 = values[N_total - 3]; f1 = values[N_total - 2]; f2 = values[N_total - 1];
    }

    // Проверка на совпадающие узлы (знаменатель не должен быть нулем)
    if (std::abs(x0 - x1) < EPS_SMALL || std::abs(x0 - x2) < EPS_SMALL || std::abs(x1 - x2) < EPS_SMALL) {
        // Совпадающие узлы, невозможно вычислить. Возвращаем 0 (как для естественного сплайна)
        return 0.0;
    }
    
    // Формула для второй производной параболы, проходящей через (x0,f0), (x1,f1), (x2,f2):
    // M = 2 * ( f0/((x0-x1)(x0-x2)) + f1/((x1-x0)(x1-x2)) + f2/((x2-x0)(x2-x1)) )
    double term0 = f0 / ((x0 - x1) * (x0 - x2));
    double term1 = f1 / ((x1 - x0) * (x1 - x2));
    double term2 = f2 / ((x2 - x0) * (x2 - x1));
    
    return 2.0 * (term0 + term1 + term2);
}


// Вычисляет 1D кубический сплайн и возвращает ПЕРВЫЕ производные в узлах.
// Граничные условия (вторые производные M0, Mn-1) вычисляются внутри экстраполяцией.
void Scene3D::compute_1D_Spline_derivatives(const double* input_vals, const double* coords, int N,
                                            double* output_first_derivatives)
{
    if (N < 2) {
        if (N == 1 && output_first_derivatives) output_first_derivatives[0] = 0.0;
        return;
    }
    if (N == 2) { // Один отрезок, производные на концах равны наклону отрезка
        double h0 = coords[1] - coords[0];
        if (std::abs(h0) < EPS_SMALL) { // Совпадающие узлы
            output_first_derivatives[0] = 0.0;
            output_first_derivatives[1] = 0.0;
        } else {
            double slope = (input_vals[1] - input_vals[0]) / h0;
            output_first_derivatives[0] = slope;
            output_first_derivatives[1] = slope;
        }
        return;
    }

    // N >= 3
    std::vector<double> h(N - 1);
    for (int i = 0; i < N - 1; ++i) {
        h[i] = coords[i+1] - coords[i];
        if (std::abs(h[i]) < EPS_SMALL) {
             // Обработка случая очень малого шага, если необходимо.
             // Для равномерной сетки и N>=3 это маловероятно, если b-a > 0.
             // Если все же произойдет, можно установить производные в 0 или выдать ошибку.
            for(int k=0; k<N; ++k) output_first_derivatives[k] = 0.0;
            std::cerr << "Warning: Very small step h[" << i << "] in 1D spline computation." << std::endl;
            return;
        }
    }

    std::vector<double> M(N); // Вторые производные M_i = P''_i(x_i)

    // 1. Граничные условия для M_0 и M_{N-1} методом экстраполяции
    M[0]   = extrapolate_M_boundary(coords, input_vals, N, true);
    M[N-1] = extrapolate_M_boundary(coords, input_vals, N, false);

    if (N == 3) { // Только M_1 неизвестна. Система из 1 уравнения.
        // h_0*M_0 + 2*(h_0+h_1)*M_1 + h_1*M_2 = 6 * ( (f_2-f_1)/h_1 - (f_1-f_0)/h_0 )
        double rhs_val = 6.0 * ( (input_vals[2]-input_vals[1])/h[1] - (input_vals[1]-input_vals[0])/h[0] );
        double diag_M1 = 2.0 * (h[0] + h[1]);
        if (std::abs(diag_M1) < EPS_SMALL) {
            M[1] = 0.0; // или другая обработка
        } else {
            M[1] = (rhs_val - h[0]*M[0] - h[1]*M[N-1]) / diag_M1; // M[N-1] здесь M[2]
        }
    } else { // N > 3, решаем систему для M_1, ..., M_{N-2}
        int system_size = N - 2; // Количество неизвестных M_1 ... M_{N-2}
        std::vector<double> tril_a_vec(system_size);    // поддиагональ (индекс от 0 до N-3, т.е. a_1 .. a_{N-2})
        std::vector<double> tril_d_vec(system_size);    // главная диагональ (d_1 .. d_{N-2})
        std::vector<double> tril_c_vec(system_size);    // наддиагональ (c_1 .. c_{N-2})
        std::vector<double> tril_b_vec(system_size);    // правая часть
        std::vector<double> M_internal_sol(system_size);// решение M_1, ..., M_{N-2}

        // Формирование системы: h_{i-1}*M_{i-1} + 2(h_{i-1}+h_i)*M_i + h_i*M_{i+1} = RHS_i
        // для i = 1, ..., N-2.
        // В нашей системе tril_X_vec[k] соответствует уравнению для M_{k+1}
        for (int k = 0; k < system_size; ++k) { // k от 0 до N-3
            int i = k + 1; // Глобальный индекс M_i (M_1, ..., M_{N-2})
            
            tril_d_vec[k] = 2.0 * (h[i-1] + h[i]);
            if (k > 0) tril_a_vec[k] = h[i-1];          // Для k=0 (M1), a[0] не используется в solveTridiagonalSystem
            if (k < system_size - 1) tril_c_vec[k] = h[i]; // Для k=N-3 (M_{N-2}), c[N-3] не используется
            
            tril_b_vec[k] = 6.0 * ( (input_vals[i+1] - input_vals[i]) / h[i] - (input_vals[i] - input_vals[i-1]) / h[i-1] );
            
            if (k == 0) { // Первое уравнение (для M_1)
                tril_b_vec[k] -= h[i-1] * M[0]; // h[0]*M[0]
            }
            if (k == system_size - 1) { // Последнее уравнение (для M_{N-2})
                tril_b_vec[k] -= h[i] * M[N-1]; // h[N-2]*M[N-1]
            }
        }
        
        // Адаптация для solveTridiagonalSystem: a[k] для (k+1)-го уравнения, d[k] для k-го, c[k] для k-го.
        // solveTridiagonalSystem ожидает: a - поддиагональ (a[k] для (k+1)-го уравнения, ее первый элемент a[0] соответствует a_1 в уравнении для x_1),
        // d - главная, c - наддиагональ.
        // В нашей tril_a_vec[k] - это коэффициент при M_k в уравнении для M_{k+1}.
        // Для solveTridiagonalSystem: a_arg[k] - коэфф. при x_k в (k+1)-м уравнении (индекс k от 0 до system_size-2)
        std::vector<double> a_arg(system_size > 1 ? system_size - 1 : 0);
        std::vector<double> c_arg(system_size > 1 ? system_size - 1 : 0);
        for(int k=0; k<system_size-1; ++k) {
            a_arg[k] = tril_a_vec[k+1]; // a_arg[0] = tril_a_vec[1] (коэфф M_1 в уравнении для M_2)
            c_arg[k] = tril_c_vec[k];   // c_arg[0] = tril_c_vec[0] (коэфф M_2 в уравнении для M_1)
        }

        // Передаем указатели на данные векторов. solveTridiagonalSystem МОЖЕТ модифицировать d и b.
        std::vector<double> d_copy = tril_d_vec;
        std::vector<double> b_copy = tril_b_vec;

        bool solved = solveTridiagonalSystem(
            system_size > 1 ? a_arg.data() : nullptr, 
            d_copy.data(), 
            system_size > 1 ? c_arg.data() : nullptr, 
            b_copy.data(), 
            M_internal_sol.data(), 
            system_size
        );

        if (solved) {
            for(int k=0; k < system_size; ++k) M[k+1] = M_internal_sol[k];
        } else {
            std::cerr << "Warning: Tridiagonal system solution failed for 1D spline M_i." << std::endl;
            for(int k=0; k < system_size; ++k) M[k+1] = 0.0; // Fallback
        }
    }

    // 3. Вычисляем первые производные d_i = P'_i(x_i) используя M_i
    // d_i = (f_{i+1}-f_i)/h_i - h_i/6 * (2*M_i + M_{i+1}) для i = 0, ..., N-2
    for (int i = 0; i < N - 1; ++i) {
        output_first_derivatives[i] = (input_vals[i+1] - input_vals[i]) / h[i] - h[i] / 6.0 * (2.0 * M[i] + M[i+1]);
    }
    // d_{N-1} = (f_{N-1}-f_{N-2})/h_{N-2} + h_{N-2}/6 * (M_{N-2} + 2*M_{N-1})
    // h[N-2] это последний шаг h_{N-2} = x_{N-1} - x_{N-2}
    output_first_derivatives[N-1] = (input_vals[N-1] - input_vals[N-2]) / h[N-2] + h[N-2] / 6.0 * (M[N-2] + 2.0 * M[N-1]);
}


void Scene3D::computeSplineData()
{
    if (nx_param < MIN_SPLINE_NODES || ny_param < MIN_SPLINE_NODES) {
        std::cerr << "Error: nx and ny must be at least " << MIN_SPLINE_NODES << " for spline interpolation." << std::endl;
        cleanupSplineData();
        return;
    }
    allocateSplineData(nx_param, ny_param);

    // 1. Инициализация узлов (равномерная сетка)
    for (int i = 0; i < current_nx; ++i) {
        spline_x_nodes[i] = a + (b - a) * i / (current_nx - 1.0);
    }
    for (int j = 0; j < current_ny; ++j) {
        spline_y_nodes[j] = c + (d - c) * j / (current_ny - 1.0);
    }

    // 2. Заполнение f_val (значения функции в узлах)
    for (int i = 0; i < current_nx; ++i) {
        for (int j = 0; j < current_ny; ++j) {
            f_val[i][j] = f(spline_x_nodes[i], spline_y_nodes[j]);
        }
    }
    if (p != 0 && current_nx > 0 && current_ny > 0) {
        f_val[current_nx/2][current_ny/2] += p * 0.1 * fmax_orig;
    }

    double* temp_values = new double[std::max(current_nx, current_ny)];
    double* temp_derivs = new double[std::max(current_nx, current_ny)];

    // 3. Вычисление fx_val (df/dx в узлах)
    for (int j = 0; j < current_ny; ++j) { // Для каждой строки Y
        for (int i = 0; i < current_nx; ++i) {
            temp_values[i] = f_val[i][j];
        }
        compute_1D_Spline_derivatives(temp_values, spline_x_nodes, current_nx, temp_derivs);
        for (int i = 0; i < current_nx; ++i) {
            fx_val[i][j] = temp_derivs[i];
        }
    }

    // 4. Вычисление fy_val (df/dy в узлах)
    for (int i = 0; i < current_nx; ++i) { // Для каждого столбца X
        for (int j = 0; j < current_ny; ++j) {
            temp_values[j] = f_val[i][j];
        }
        compute_1D_Spline_derivatives(temp_values, spline_y_nodes, current_ny, temp_derivs);
        for (int j = 0; j < current_ny; ++j) {
            fy_val[i][j] = temp_derivs[j];
        }
    }

    // 5. Вычисление fxy_val (d2f/dxdy, дифференцируем fx_val по Y)
    for (int i = 0; i < current_nx; ++i) { // Для каждого столбца X
        for (int j = 0; j < current_ny; ++j) {
            temp_values[j] = fx_val[i][j]; // fx_val теперь входные данные
        }
        compute_1D_Spline_derivatives(temp_values, spline_y_nodes, current_ny, temp_derivs);
        for (int j = 0; j < current_ny; ++j) {
            fxy_val[i][j] = temp_derivs[j];
        }
    }

    delete[] temp_values;
    delete[] temp_derivs;
}

double Scene3D::evaluateSpline(double x, double y)
{
    if (!f_val || current_nx < 2 || current_ny < 2) return 0.0; // Need at least one cell

    int cell_i = -1, cell_j = -1;

    if (x <= spline_x_nodes[0]) cell_i = 0;
    else if (x >= spline_x_nodes[current_nx-1]) cell_i = current_nx - 2;
    else {
        for (int k = 0; k < current_nx - 1; ++k) {
            if (x >= spline_x_nodes[k] && x <= spline_x_nodes[k+1]) {
                cell_i = k; break;
            }
        }
    }
    if (cell_i == -1) cell_i = 0; 

    if (y <= spline_y_nodes[0]) cell_j = 0;
    else if (y >= spline_y_nodes[current_ny-1]) cell_j = current_ny - 2;
    else {
        for (int k = 0; k < current_ny - 1; ++k) {
            if (y >= spline_y_nodes[k] && y <= spline_y_nodes[k+1]) {
                cell_j = k; break;
            }
        }
    }
    if (cell_j == -1) cell_j = 0;

    int i0 = cell_i;     int i1 = cell_i + 1;
    int j0 = cell_j;     int j1 = cell_j + 1;

    double hx = spline_x_nodes[i1] - spline_x_nodes[i0];
    double hy = spline_y_nodes[j1] - spline_y_nodes[j0];

    if (fabs(hx) < EPS_SMALL) hx = (hx > 0 ? EPS_SMALL : -EPS_SMALL); // Avoid division by zero, preserve sign if possible
    if (fabs(hy) < EPS_SMALL) hy = (hy > 0 ? EPS_SMALL : -EPS_SMALL);


    double s = (x - spline_x_nodes[i0]) / hx;
    double t = (y - spline_y_nodes[j0]) / hy;
    
    s = std::max(0.0, std::min(1.0, s)); // Clamp s and t to [0,1] for robustness
    t = std::max(0.0, std::min(1.0, t));


    double G[4][4];
    G[0][0] = f_val[i0][j0];        G[0][1] = f_val[i0][j1];        G[0][2] = fy_val[i0][j0] * hy;    G[0][3] = fy_val[i0][j1] * hy;
    G[1][0] = f_val[i1][j0];        G[1][1] = f_val[i1][j1];        G[1][2] = fy_val[i1][j0] * hy;    G[1][3] = fy_val[i1][j1] * hy;
    G[2][0] = fx_val[i0][j0] * hx;   G[2][1] = fx_val[i0][j1] * hx;   G[2][2] = fxy_val[i0][j0] * hx * hy; G[2][3] = fxy_val[i0][j1] * hx * hy;
    G[3][0] = fx_val[i1][j0] * hx;   G[3][1] = fx_val[i1][j1] * hx;   G[3][2] = fxy_val[i1][j0] * hx * hy; G[3][3] = fxy_val[i1][j1] * hx * hy;

    double C_h[4][4] = {
        { 1,  0,  0,  0}, { 0,  0,  1,  0},
        {-3,  3, -2, -1}, { 2, -2,  1,  1}
    };
    double C_h_T[4][4];
    for(int r=0; r<4; ++r) for(int col=0; col<4; ++col) C_h_T[r][col] = C_h[col][r];

    double Temp_Matrix[4][4] = {{0}};
    for(int r=0; r<4; ++r) for(int c_p=0; c_p<4; ++c_p) for(int k=0; k<4; ++k) Temp_Matrix[r][c_p] += C_h[r][k] * G[k][c_p];
    
    double Coeff_Matrix[4][4] = {{0}};
    for(int r=0; r<4; ++r) for(int c_p=0; c_p<4; ++c_p) for(int k=0; k<4; ++k) Coeff_Matrix[r][c_p] += Temp_Matrix[r][k] * C_h_T[k][c_p];

    double S_vec[4] = {1, s, s*s, s*s*s};
    double T_vec[4] = {1, t, t*t, t*t*t};
    double Temp_S_Coeff[4] = {0};
    for(int c_p=0; c_p<4; ++c_p) for(int k=0; k<4; ++k) Temp_S_Coeff[c_p] += S_vec[k] * Coeff_Matrix[k][c_p];
    
    double result = 0;
    for(int k=0; k<4; ++k) result += Temp_S_Coeff[k] * T_vec[k];
    return result;
}

double Scene3D::errorSpline(double x, double y) {
    if (!f_val || !f) return 0.0;
    return fabs(evaluateSpline(x, y) - f(x, y));
}

void Scene3D::initializeGL() {
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
}

void Scene3D::resizeGL(int nWidth, int nHeight) {
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;
    if (nWidth >= nHeight) glOrtho(-1.0 / ratio, 1.0 / ratio, -1.0, 1.0, -10.0, 1.0);
    else glOrtho(-1.0, 1.0, -1.0 * ratio, 1.0 * ratio, -10.0, 1.0);
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

void Scene3D::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();
    glScalef(nSca, nSca, nSca);
    glTranslatef(0.0f, zTra, 0.0f);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);
    glRotatef(zRot, 0.0f, 0.0f, 1.0f);

    char str_buf[100];
    glColor3d(0,0,1);
    renderText(0, 20, f_name);
    if (fmax_orig > 1e5) std::snprintf(str_buf, 100, "Max|f_orig| = %.2e", fmax_orig);
    else std::snprintf(str_buf, 100, "Max|f_orig| = %.4f", fmax_orig);
    renderText(0, 35, str_buf);
    std::snprintf(str_buf, 100, "nx, ny ~ %d", (current_nx + current_ny)/2);
    renderText(0, 80, str_buf);
    std::snprintf(str_buf, 100, "angle = %d", repeat);
    renderText(0, 65, str_buf);
    std::snprintf(str_buf, 100, "error p = %d", p);
    renderText(0, 50, str_buf);

    drawAxis();
    if (gr == 0) drawOriginalFunction();
    else if (gr == 1) {
        drawSplineApproximation(); // Calls evaluateSpline, updates max_metod
        if (max_metod > 1e5) std::snprintf(str_buf, 100, "Max|f_spline| = %.2e", max_metod);
        else std::snprintf(str_buf, 100, "Max|f_spline| = %.4f", max_metod);
        renderText(0, 95, str_buf);
    } else if (gr == 2) drawErrorFunction();
}

void Scene3D::mousePressEvent(QMouseEvent *pe) { ptrMousePosition = pe->pos(); }
void Scene3D::mouseReleaseEvent(QMouseEvent *pe) { UNUSED(pe); }
void Scene3D::mouseMoveEvent(QMouseEvent *pe) {
    xRot += 180 / nSca * (GLfloat)(pe->y() - ptrMousePosition.y()) / height();
    zRot += 180 / nSca * (GLfloat)(pe->x() - ptrMousePosition.x()) / width();
    ptrMousePosition = pe->pos(); updateGL();
}
void Scene3D::wheelEvent(QWheelEvent *pe) {
    if ((pe->delta()) > 0) scale_plus(); else if ((pe->delta()) < 0) scale_minus();
    updateGL();
}
void Scene3D::keyPressEvent(QKeyEvent *pe) {
    switch (pe->key()) {
        case Qt::Key_Escape: close(); break;
        case Qt::Key_Up: rotate_up(); break; case Qt::Key_Down: rotate_down(); break;
        case Qt::Key_Left: rotate_left(); break; case Qt::Key_Right: rotate_right(); break;
        case Qt::Key_1: change_graph(); break; case Qt::Key_0: change_func(); break;
        case Qt::Key_2: increase_param_scale(); break; case Qt::Key_3: decrease_param_scale(); break;
        case Qt::Key_4: increase_nx_ny_param(); break; case Qt::Key_5: decrease_nx_ny_param(); break;
        case Qt::Key_6: increase_error_p(); break; case Qt::Key_7: decrease_error_p(); break;
        case Qt::Key_8: plus_angle(); break; case Qt::Key_9: minus_angle(); break;
        case Qt::Key_Space: defaultScene(); break;
    }
    updateGL();
}

void Scene3D::scale_plus() { nSca *= 1.1; } void Scene3D::scale_minus() { nSca /= 1.1; }
void Scene3D::rotate_up() { xRot += 1.0; } void Scene3D::rotate_down() { xRot -= 1.0; }
void Scene3D::rotate_left() { zRot += 1.0; } void Scene3D::rotate_right() { zRot -= 1.0; }
void Scene3D::p_angle() { zRot += 15; repeat = (repeat+15)%360; }
void Scene3D::m_angle() { zRot -= 15; repeat = (repeat-15+360)%360; }
void Scene3D::defaultScene() { xRot = -90; yRot = 0; zRot = 0; zTra = 0; nSca = 1; }

void Scene3D::drawAxis() {
    glLineWidth(3.0f);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f); glBegin(GL_LINES); glVertex3f(-1.0f, 0.0f, 0.0f); glVertex3f(1.0f, 0.0f, 0.0f); glEnd();
    glColor4f(0.0f, 1.0f, 0.0f, 1.0f); glBegin(GL_LINES); glVertex3f(0.0f, -1.0f, 0.0f); glVertex3f(0.0f, 1.0f, 0.0f); glEnd();
    glColor4f(0.0f, 0.0f, 1.0f, 1.0f); glBegin(GL_LINES); glVertex3f(0.0f, 0.0f, -1.0f); glVertex3f(0.0f, 0.0f, 1.0f); glEnd();
}

void Scene3D::drawOriginalFunction() {
    if (!f) return;
    int N_draw = 50;
    double dx_draw = (b - a) / N_draw; double dy_draw = (d - c) / N_draw;
    glLineWidth(1.5f); qglColor(Qt::black); glBegin(GL_LINES);
    for (int i = 0; i <= N_draw; ++i) {
        double cur_x = a + i * dx_draw;
        double prev_y = c, prev_z = f(cur_x, prev_y);
        for (int j = 1; j <= N_draw; ++j) {
            double cur_y = c + j * dy_draw, cur_z = f(cur_x, cur_y);
            glVertex3d(cur_x, prev_y, prev_z); glVertex3d(cur_x, cur_y, cur_z);
            prev_y = cur_y; prev_z = cur_z;
        }
    }
    for (int j = 0; j <= N_draw; ++j) {
        double cur_y = c + j * dy_draw;
        double prev_x = a, prev_z = f(prev_x, cur_y);
        for (int i = 1; i <= N_draw; ++i) {
            double cur_x = a + i * dx_draw, cur_z = f(cur_x, cur_y);
            glVertex3d(prev_x, cur_y, prev_z); glVertex3d(cur_x, cur_y, cur_z);
            prev_x = cur_x; prev_z = cur_z;
        }
    }
    glEnd();
}

void Scene3D::drawSplineApproximation() {
    if (!f_val) return;
    int N_draw = 50;
    double dx_draw = (b - a) / N_draw; double dy_draw = (d - c) / N_draw;
    glLineWidth(1.0f); qglColor(Qt::cyan); glBegin(GL_LINES);
    max_metod = 0;
    for (int i = 0; i <= N_draw; ++i) {
        double cur_x = a + i * dx_draw;
        double prev_y = c, prev_z = evaluateSpline(cur_x, prev_y);
        if (fabs(prev_z) > max_metod) max_metod = fabs(prev_z);
        for (int j = 1; j <= N_draw; ++j) {
            double cur_y = c + j * dy_draw, cur_z = evaluateSpline(cur_x, cur_y);
            if (fabs(cur_z) > max_metod) max_metod = fabs(cur_z);
            glVertex3d(cur_x, prev_y, prev_z); glVertex3d(cur_x, cur_y, cur_z);
            prev_y = cur_y; prev_z = cur_z;
        }
    }
    for (int j = 0; j <= N_draw; ++j) {
        double cur_y = c + j * dy_draw;
        double prev_x = a, prev_z = evaluateSpline(prev_x, cur_y);
        if (fabs(prev_z) > max_metod) max_metod = fabs(prev_z);
        for (int i = 1; i <= N_draw; ++i) {
            double cur_x = a + i * dx_draw, cur_z = evaluateSpline(cur_x, cur_y);
            if (fabs(cur_z) > max_metod) max_metod = fabs(cur_z);
            glVertex3d(prev_x, cur_y, prev_z); glVertex3d(cur_x, cur_y, cur_z);
            prev_x = cur_x; prev_z = cur_z;
        }
    }
    glEnd();
}

void Scene3D::drawErrorFunction() {
    if (!f_val || !f) return;
    int N_draw = 50;
    double dx_draw = (b - a) / N_draw; double dy_draw = (d - c) / N_draw;
    glLineWidth(1.0f); qglColor(Qt::magenta); glBegin(GL_LINES);
    // double max_err_val = 0.0; // Optional: to display max error
    for (int i = 0; i <= N_draw; ++i) {
        double cur_x = a + i * dx_draw;
        double prev_y = c, prev_z_err = errorSpline(cur_x, prev_y);
        // if (prev_z_err > max_err_val) max_err_val = prev_z_err;
        for (int j = 1; j <= N_draw; ++j) {
            double cur_y = c + j * dy_draw, cur_z_err = errorSpline(cur_x, cur_y);
            // if (cur_z_err > max_err_val) max_err_val = cur_z_err;
            glVertex3d(cur_x, prev_y, prev_z_err); glVertex3d(cur_x, cur_y, cur_z_err);
            prev_y = cur_y; prev_z_err = cur_z_err;
        }
    }
    for (int j = 0; j <= N_draw; ++j) {
        double cur_y = c + j * dy_draw;
        double prev_x = a, prev_z_err = errorSpline(prev_x, cur_y);
        // if (prev_z_err > max_err_val) max_err_val = prev_z_err;
        for (int i = 1; i <= N_draw; ++i) {
            double cur_x = a + i * dx_draw, cur_z_err = errorSpline(cur_x, cur_y);
            // if (cur_z_err > max_err_val) max_err_val = cur_z_err;
            glVertex3d(prev_x, cur_y, prev_z_err); glVertex3d(cur_x, cur_y, cur_z_err);
            prev_x = cur_x; prev_z_err = cur_z_err;
        }
    }
    glEnd();
    // char err_str[100];
    // std::snprintf(err_str, 100, "Max|Error| = %.2e", max_err_val);
    // renderText(0, 110, err_str);
}
