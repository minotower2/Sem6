#include <cstdio>
#include <limits>
#include <functional>
#include <cmath>
#include <QPainter>
#include <QPointF>
#include <algorithm>
#include <vector> 

#include "window.h"
#include "functions.h"
#include "chebyshev_approximation.h" 
#include "hermite_spline_approximation.h" 
#include "parabolic_approximation.h"

Window::Window(QWidget *parent, double a_in, double b_in, int n_in, int k_in)
    : QWidget(parent), a(a_in), show_a(a_in), b(b_in), show_b(b_in), n(n_in), k(k_in) {
    
    x.reset(new double[n]); 
    f.reset(new double[n]);
    d.reset(new double[n]);
    xi.reset(new double[n+1]);

    function = get_function(k);
    derivative = get_derivative(k);

    calculate_points();

    chebyshev_c.reset(new double[n]);
    make_chebyshev_coefficients(n, this->a, this->b, function, chebyshev_c.get());

    hermite_c.reset(new double[4 * n]); 
    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get());

    parabolic_c.reset(new double[3*n]);
    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());
    
    calculate_max_abs_f();
    calculate_min_max();
    init_status_bar();
}

void Window::calculate_points() { // Для узлов Эрмита (равномерная сетка)
    double x_0_val;
    double h_left, h_right;
    if (n < 1) return;
    if (n == 1) {
        x[0] = a;
        f[0] = function(a);
        d[0] = derivative(a);
	xi[0] = a;
        return;
    }
    for (int i = 0; i < n; i++) {
        x_0_val = a + i * (b - a) / (n - 1);
        x[i] = x_0_val;
        f[i] = function(x_0_val);
        d[i] = derivative(x_0_val);
    }
    h_left = x[1] - x[0];
    h_right = x[n-1] - x[n-2];
    for (int i = 1; i < n; i++) {
        xi[i] = 0.5*(x[i] + x[i-1]);
    }
    xi[0] = x[0] - 0.5* h_left;
    xi[n] = x[n-1] + 0.5 * h_right;
}

void Window::calculate_max_abs_f() {
    double x_0_val, abs_y_0;
    double delta_x = (b - a) / (width() > 0 ? static_cast<double>(width()) : 1000.0);
    if (delta_x <= 0 && b > a) delta_x = (b-a) / 1000.0; // Fallback if width is 0 or too small
    else if (delta_x <= 0) { // If b <= a
        max_abs_f = std::abs(function(a)); 
        return;
    }

    max_abs_f = std::fabs(function(a));
    for (x_0_val = a + delta_x; x_0_val - b < std::numeric_limits<double>::epsilon(); x_0_val += delta_x) {
        abs_y_0 = std::fabs(function(x_0_val));
        if (abs_y_0 > max_abs_f) {
            max_abs_f = abs_y_0;
        }
    }
    abs_y_0 = std::fabs(function(b));
    if (abs_y_0 > max_abs_f) {
        max_abs_f = abs_y_0;
    }
}

void Window::calculate_min_max() {
    double x_0_val, y_0_val;
    double delta_x = (show_b - show_a) / (width() > 0 ? static_cast<double>(width()) : 1000.0);
    if (delta_x <= 0 && show_b > show_a) delta_x = (show_b - show_a) / 1000.0;
    else if (delta_x <= 0) {
        min_y = function(show_a) -1.0; max_y = function(show_a) +1.0;
        max_abs_F = std::fmax(std::fabs(min_y), std::fabs(max_y));
        char str_f_loc[24];
        std::snprintf(str_f_loc, 24, "max{|F|} = %9.2e", max_abs_F);
        max_abs_F_changed(QString(str_f_loc));
        return;
    }

    bool is_discrepancy_mode = (mode >= 4 && mode <= 6);
    bool first_val_set = false;
    min_y = std::numeric_limits<double>::infinity();
    max_y = -std::numeric_limits<double>::infinity();

    // Режимы:
    // 0: Исходная функция
    // 1: Исходная + Чебышев
    // 2: Исходная + Эрмит
    // 3: Исходная + параболическая
    // 4: Исходная + Чебышев + Эрмит + параболическая
    // 5: Невязка Чебышева
    // 6: Невязка Эрмита
    // 7: Невязка параболической
    // 8: Три невязки
    
    for (x_0_val = show_a; x_0_val <= show_b + std::numeric_limits<double>::epsilon(); x_0_val += delta_x) {
        // Исходная функция (если не чисто режим невязок без исходной)
        if (mode <= 4) { // Режимы 0, 1, 2, 3, 4 включают исходную функцию
            y_0_val = function(x_0_val); 
             if (!std::isnan(y_0_val)) {
                if (!first_val_set) { min_y = max_y = y_0_val; first_val_set = true; }
                else { min_y = std::min(min_y, y_0_val); max_y = std::max(max_y, y_0_val); }
            }
        }

        // Чебышев
        if (mode == 1 || mode == 4 || mode == 5 || mode == 8) {
            if (n >= 1 && chebyshev_c && !std::isnan(chebyshev_c[0])) {
                y_0_val = (mode == 1 || mode == 4) ? calculate_chebyshev_approximation(x_0_val, n, this->a, this->b, chebyshev_c.get())
                                               : calculate_chebyshev_discrepancy(x_0_val, function, n, this->a, this->b, chebyshev_c.get());
                if (!std::isnan(y_0_val)) {
                    if (!first_val_set || (is_discrepancy_mode && mode > 4) ) { min_y = max_y = y_0_val; first_val_set = true; }
                    else { min_y = std::min(min_y, y_0_val); max_y = std::max(max_y, y_0_val); }
                }
            }
        }
        // Эрмит
        if (mode == 2 || mode == 4 || mode == 6 || mode == 8) {
            if (n >= 2 && hermite_c && !std::isnan(hermite_c[0])) {
                y_0_val = (mode == 2 || mode == 4) ? calculate_cubic_hermite_approximation(x_0_val, n, x.get(), hermite_c.get())
                                               : calculate_cubic_hermite_discrepancy(x_0_val, function, n, x.get(), hermite_c.get());
                if (!std::isnan(y_0_val)) {
                    if (!first_val_set || (is_discrepancy_mode && mode > 4) ) { min_y = max_y = y_0_val; first_val_set = true; }
                    else { min_y = std::min(min_y, y_0_val); max_y = std::max(max_y, y_0_val); }
                }
            }
        }
        // Параболический
        if (mode == 3 || mode == 4 || mode == 7 || mode == 8) {
            if (n >= 2 && hermite_c && !std::isnan(hermite_c[0])) {
                y_0_val = (mode == 3 || mode == 4) ? calculate_parabolic_spline_approximation(x_0_val, n, x.get(), xi.get(), parabolic_c.get())
                                               : calculate_parabolic_spline_discrepancy(x_0_val, function, n, x.get(), xi.get(), hermite_c.get());
                if (!std::isnan(y_0_val)) {
                    if (!first_val_set || (is_discrepancy_mode && mode > 4) ) { min_y = max_y = y_0_val; first_val_set = true; }
                    else { min_y = std::min(min_y, y_0_val); max_y = std::max(max_y, y_0_val); }
                }
            }
        }
    }
    if (!first_val_set) { 
        min_y = -1.0; max_y = 1.0;
    }

    max_abs_F = std::fmax(std::fabs(min_y), std::fabs(max_y));
    char str_f_loc[24];
    std::snprintf(str_f_loc, 24, "max{|F|} = %9.2e", max_abs_F);
    max_abs_F_changed(QString(str_f_loc));

    double delta_y_padding;
    if (max_y - min_y < std::numeric_limits<double>::epsilon() ) { 
        delta_y_padding = 0.1; 
    } else {
        delta_y_padding = 0.05 * (max_y - min_y);
    }
     if (!std::isinf(min_y) && !std::isinf(max_y) && !std::isnan(min_y) && !std::isnan(max_y)) {
        min_y -= delta_y_padding;
        max_y += delta_y_padding;
     } else { 
        min_y = -1.0; max_y = 1.0; 
     }
}


void Window::init_status_bar() {
    function_description_changed(QString(get_function_description(k)));
    char str_buf[24];
    std::snprintf(str_buf, 16, "n = %d", n);
    number_of_points_changed(QString(str_buf));
    std::snprintf(str_buf, 16, "s = %d", s);
    scale_changed(QString(str_buf));
    std::snprintf(str_buf, 16, "p = %d", p);
    distortion_changed(QString(str_buf));
    std::snprintf(str_buf, 24, "max{|F|} = %9.2e", max_abs_F);
    max_abs_F_changed(QString(str_buf));
}

void Window::change_function() {
    k++;
    if (k > 6) k = 0;
    function = get_function(k);
    derivative = get_derivative(k);

    
    calculate_points(); // x, f обновлены оригинальными значениями

    p = 0; 
    
    chebyshev_c.reset(new double[n]); 
    make_chebyshev_coefficients(n, this->a, this->b, function, chebyshev_c.get()); 

    hermite_c.reset(new double[4*n]); 
    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get());

    parabolic_c.reset(new double[3*n]);
    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());
    
    calculate_max_abs_f(); 
    calculate_min_max();
    function_description_changed(QString(get_function_description(k)));
    
    char str_p_loc[16]; 
    std::snprintf(str_p_loc, 16, "p = %d", p);
    distortion_changed(QString(str_p_loc));
    
    update();
}

void Window::change_graph() {
    mode++;
    if (mode >= total_modes) { 
        mode = 0;
    }
    calculate_min_max(); 
    update();
}

void Window::increase_scale() {
    s++;
    double center = (show_a + show_b) * 0.5;
    show_a = (show_a + center) * 0.5;
    show_b = (show_b + center) * 0.5;
    char str_s_loc[16];
    std::snprintf(str_s_loc, 16, "s = %d", s);
    scale_changed(QString(str_s_loc));
    calculate_min_max();
    update();
}

void Window::decrease_scale() {
    s--;
    double center = (show_a + show_b) * 0.5;
    show_a = 2 * show_a - center;
    show_b = 2 * show_b - center;
    char str_s_loc[16];
    std::snprintf(str_s_loc, 16, "s = %d", s);
    scale_changed(QString(str_s_loc));
    calculate_min_max();
    update();
}

void Window::increase_points() {
    n *= 2;
    if (n <= 0) n = 1;

    x.reset(new double[n]);
    f.reset(new double[n]);
    d.reset(new double[n]);
    xi.reset(new double[n+1]);
    calculate_points(); 

    if (p != 0 && n >= 1) { 
        f[n / 2] += static_cast<double>(p) * 0.1 * max_abs_f;
    }
    
    chebyshev_c.reset(new double[n]);
    make_chebyshev_coefficients(n, this->a, this->b, function, chebyshev_c.get()); 

    hermite_c.reset(new double[4*n]);
    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get());

    parabolic_c.reset(new double[3*n]);
    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());

    calculate_min_max();
    char str_n_loc[16];
    std::snprintf(str_n_loc, 16, "n = %d", n);
    number_of_points_changed(QString(str_n_loc));
    update();
}

void Window::decrease_points() {
    n /= 2;
    if (n <= 0) n = 1;

    x.reset(new double[n]);
    f.reset(new double[n]);
    d.reset(new double[n]);
    xi.reset(new double[n+1]);
    calculate_points();

    if (p != 0 && n >= 1) {
        f[n / 2] += static_cast<double>(p) * 0.1 * max_abs_f;
    }

    chebyshev_c.reset(new double[n]);
    make_chebyshev_coefficients(n, this->a, this->b, function, chebyshev_c.get());

    hermite_c.reset(new double[4*n]);
    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get());


    parabolic_c.reset(new double[3*n]);
    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());

    calculate_min_max();
    char str_n_loc[16];
    std::snprintf(str_n_loc, 16, "n = %d", n);
    number_of_points_changed(QString(str_n_loc));
    update();
}

void Window::add_distortion() {
    if (n < 1) return;
    p++;
    f[n / 2] += 0.1 * max_abs_f; 

    // Пересчитываем Бесселя, т.к. он зависит от f.get()
    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get()); 
    // Чебышев не пересчитывается, он зависит от function.

    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());

    char str_p_dist[16];
    std::snprintf(str_p_dist, 16, "p = %d", p);
    distortion_changed(QString(str_p_dist));
    calculate_min_max();
    update();
}

void Window::subtract_distortion() {
    if (n < 1) return;
    p--;
    f[n / 2] -= 0.1 * max_abs_f;

    make_cubic_hermite_coefficients(n, x.get(), f.get(), d.get(), hermite_c.get());

    make_parabolic_spline_coefficients(n, x.get(), f.get(), xi.get(), parabolic_c.get());

    char str_p_dist[16];
    std::snprintf(str_p_dist, 16, "p = %d", p);
    distortion_changed(QString(str_p_dist));
    calculate_min_max();
    update();
}

QPointF Window::l2g(double x_coord, double y_coord) {
    double x_gl, y_gl;
    if (std::abs(show_b - show_a) < std::numeric_limits<double>::epsilon()) {
        x_gl = width() / 2.0;
    } else {
        x_gl = (x_coord - show_a) / (show_b - show_a) * width();
    }

    if (std::abs(max_y - min_y) < std::numeric_limits<double>::epsilon() || 
        std::isinf(min_y) || std::isinf(max_y) || std::isnan(min_y) || std::isnan(max_y)) {
        if (fabs(y_coord) < 1e-12) y_gl = height() / 2.0;
        else if (y_coord > 0) y_gl = height() / 4.0; 
        else y_gl = 3 * height() / 4.0;          
    } else {
        y_gl = (max_y - y_coord) / (max_y - min_y) * height();
    }
    return QPointF(x_gl, y_gl);
}

template<typename Function>
void Window::drawGraph(QPainter& painter, Function graph_func) {
    double delta_x_draw = (show_b - show_a) / width();
     if (delta_x_draw <= 0 && show_b > show_a) delta_x_draw = (show_b - show_a) / 1000.0;
     else if (delta_x_draw <= 0) return; 

    double x1_draw, y1_draw, x2_draw, y2_draw;
    x1_draw = show_a;
    y1_draw = graph_func(x1_draw);

    for (x2_draw = show_a + delta_x_draw; x2_draw <= show_b + std::numeric_limits<double>::epsilon(); x2_draw += delta_x_draw) {
        y2_draw = graph_func(x2_draw);
        if (!std::isnan(y1_draw) && !std::isnan(y2_draw) &&
            !std::isinf(y1_draw) && !std::isinf(y2_draw)) { 
             painter.drawLine(l2g(x1_draw, y1_draw), l2g(x2_draw, y2_draw));
        }
        x1_draw = x2_draw;
        y1_draw = y2_draw;
    }
    
    x2_draw = show_b;
    y2_draw = graph_func(x2_draw);
    if (!std::isnan(y1_draw) && !std::isnan(y2_draw) &&
        !std::isinf(y1_draw) && !std::isinf(y2_draw)) {
        painter.drawLine(l2g(x1_draw, y1_draw), l2g(x2_draw, y2_draw));
    }
}

void Window::paintEvent(QPaintEvent* event_paint) {
    (void)event_paint;
    QPainter painter(this);
    QPen pen_black(Qt::black, 2);
    QPen pen_red(QColor("red"), 3);         // Исходная функция
    // pen_blue удален (Ньютон)
    // pen_green удален (старый сплайн)
    QPen pen_magenta(QColor("magenta"), 3);  // Чебышев
    QPen pen_cyan(QColor("darkCyan"), 3);     // Эрмит

    char str_axis1[16], str_axis2[16];

    // Режимы:
    // 0: Исходная функция
    // 1: Исходная + Чебышев
    // 2: Исходная + Бессель
    // 3: Исходная + Чебышев + Бессель
    // 4: Невязка Чебышева
    // 5: Невязка Эрмита
    // 6: Обе невязки

    // Режимы:
    // 0: Исходная функция
    // 1: Исходная + Чебышев
    // 2: Исходная + Эрмит
    // 3: Исходная + параболическая
    // 4: Исходная + Чебышев + Эрмит + параболическая
    // 5: Невязка Чебышева
    // 6: Невязка Эрмита
    // 7: Невязка параболической
    // 8: Три невязки

    // Отрисовка исходной функции
    if (mode <= 3) { // Режимы 0, 1, 2, 3
        painter.setPen(pen_red);
        drawGraph(painter, function);
    }
    
    // Чебышев
    if (mode == 1 || mode == 4) { // Аппроксимация
        if ((n >= 1 && chebyshev_c && !std::isnan(chebyshev_c[0])) && (n <= 100)) {
            painter.setPen(pen_magenta);
            drawGraph(painter, std::bind(calculate_chebyshev_approximation, std::placeholders::_1, n, this->a, this->b, chebyshev_c.get()));
        }
    } else if (mode == 5 || mode == 8) { // Невязка
        if (n >= 1 && chebyshev_c && !std::isnan(chebyshev_c[0])) {
            painter.setPen(pen_magenta);
            drawGraph(painter, std::bind(calculate_chebyshev_discrepancy, std::placeholders::_1, function, n, this->a, this->b, chebyshev_c.get()));
        }
    }

    // Эрмит сплайн
    if (mode == 2 || mode == 4) { // Аппроксимация
        if (n >= 2 && hermite_c && !std::isnan(hermite_c[0]) ) { // Эрмит требует n>=2
            painter.setPen(pen_cyan);
            drawGraph(painter, std::bind(calculate_cubic_hermite_approximation, std::placeholders::_1, n, x.get(), hermite_c.get()));
        }
    } else if (mode == 6 || mode == 8) { // Невязка
         if (n >= 2 && hermite_c && !std::isnan(hermite_c[0]) ) {
            painter.setPen(pen_cyan);
            drawGraph(painter, std::bind(calculate_cubic_hermite_discrepancy, std::placeholders::_1, function, n, x.get(), hermite_c.get()));
        }
    }

    // Параболический сплайн
    if (mode == 3 || mode == 4) { // Аппроксимация
        if (n >= 2 && parabolic_c && !std::isnan(parabolic_c[0]) ) {
            painter.setPen(pen_black);
            drawGraph(painter, std::bind(calculate_parabolic_spline_approximation, std::placeholders::_1, n, x.get(), xi.get(), parabolic_c.get()));
        }
    } else if (mode == 7 || mode == 8) { // Невязка
         if (n >= 2 && parabolic_c && !std::isnan(parabolic_c[0]) ) {
            painter.setPen(pen_black);
            drawGraph(painter, std::bind(calculate_parabolic_spline_discrepancy, std::placeholders::_1, function, n, x.get(), xi.get(), parabolic_c.get()));
        }
    }
    painter.setPen(pen_black);
    painter.drawLine(l2g(show_a, 0), l2g(show_b, 0)); 
    painter.drawLine(l2g(0, min_y), l2g(0, max_y)); 

    double label_offset_x_factor = 0.01;
    double label_offset_y_factor = 0.03;
    double tick_len_y_factor = 0.02;
    double tick_len_x_factor = 0.02;

    double show_a_lbl = show_a + label_offset_x_factor * (show_b - show_a);
    double show_b_lbl = show_b - label_offset_x_factor * (show_b - show_a);
    double min_y_lbl = min_y + label_offset_x_factor * (max_y - min_y); 
    double max_y_lbl = max_y - label_offset_x_factor * (max_y - min_y);

    painter.drawLine(l2g(show_a_lbl, -tick_len_y_factor * (max_y - min_y)), l2g(show_a_lbl, tick_len_y_factor * (max_y - min_y)));
    painter.drawLine(l2g(show_b_lbl, -tick_len_y_factor * (max_y - min_y)), l2g(show_b_lbl, tick_len_y_factor * (max_y - min_y)));
    painter.drawLine(l2g(-tick_len_x_factor * (show_b - show_a), min_y_lbl), l2g(tick_len_x_factor * (show_b - show_a), min_y_lbl));
    painter.drawLine(l2g(-tick_len_x_factor * (show_b - show_a), max_y_lbl), l2g(tick_len_x_factor * (show_b - show_a), max_y_lbl));

    bool use_exp_x = (std::abs(show_b - show_a) < 1e-2 || std::abs(show_b - show_a) > 1e3 || std::abs(show_a_lbl) > 1e3 || std::abs(show_b_lbl) > 1e3);
    std::snprintf(str_axis1, 16, use_exp_x ? "%9.2e" : "%.2f", show_a_lbl);
    std::snprintf(str_axis2, 16, use_exp_x ? "%9.2e" : "%.2f", show_b_lbl);
    
    QPointF p_x1 = l2g(show_a_lbl, -label_offset_y_factor * (max_y - min_y)); 
    QPointF p_x2 = l2g(show_b_lbl, -label_offset_y_factor * (max_y - min_y));
    painter.drawText(p_x1, str_axis1);
    painter.drawText(p_x2 - QPointF(painter.fontMetrics().horizontalAdvance(str_axis2), 0), str_axis2);

    bool use_exp_y = (std::abs(max_y - min_y) < 1e-2 || std::abs(max_y - min_y) > 1e3 || std::abs(min_y_lbl) > 1e3 || std::abs(max_y_lbl) > 1e3);
    std::snprintf(str_axis1, 16, use_exp_y ? "%9.2e" : "%.2f", min_y_lbl); 
    std::snprintf(str_axis2, 16, use_exp_y ? "%9.2e" : "%.2f", max_y_lbl); 

    QPointF p_y1 = l2g(-label_offset_y_factor * (show_b - show_a), min_y_lbl); 
    QPointF p_y2 = l2g(-label_offset_y_factor * (show_b - show_a), max_y_lbl);
    painter.drawText(p_y1 - QPointF(painter.fontMetrics().horizontalAdvance(str_axis1), -painter.fontMetrics().ascent()/2.0), str_axis1);
    painter.drawText(p_y2 - QPointF(painter.fontMetrics().horizontalAdvance(str_axis2), -painter.fontMetrics().ascent()/2.0), str_axis2);
}
