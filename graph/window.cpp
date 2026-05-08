#include <cstdio>
#include <limits>
#include <functional>
#include <cmath>
#include <QPainter>
#include "window.h"
#include "functions.h"
#include "newton_approximation.h"
#include "spline_approximation.h"

Window::Window(QWidget *parent, double a, double b, int n, int k) : QWidget(parent), a(a), show_a(a), b(b), show_b(b), n(n), k(k) {
	x.reset(new double[n]);
	f.reset(new double[n]);
	d.reset(new double[n]);
	function = get_function(k);
	derivative = get_derivative(k);
	calculate_points();
	if(n <= 50){
		newton_c.reset(new double[2 * n]);
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	}
	spline_c.reset(new double[4 * n]);
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_max_abs_f();
	calculate_min_max();
}

void Window::calculate_points(){
	double x_0;
	if(n == 1){
		x[0] = a;
		f[0] = function(a);
		d[0] = derivative(a);
		return;
	}
	for(int i = 0;i < n;i++){
		x_0 = a + i * (b - a) / (n - 1);
		x[i] = x_0;
		f[i] = function(x_0);
		d[i] = derivative(x_0);
	}
}

void Window::calculate_max_abs_f(){
	double x_0, abs_y_0;
	double delta_x = (b - a) / width();
	max_abs_f = std::fabs(function(a));
	for(x_0 = a + delta_x; x_0 - b < std::numeric_limits<double>::epsilon(); x_0 += delta_x){
		abs_y_0 = std::fabs(function(x_0));
		if(abs_y_0 > max_abs_f){
			max_abs_f = abs_y_0;
		}
	}
}

// calculate min and max for current function
void Window::calculate_min_max(){
	double x_0, y_0;
	double delta_x = (show_b - show_a) / width();
	if(mode != 3){
		max_y = min_y = function(show_a);
	} else {
		max_y = min_y = calculate_spline_discrepancy(show_a, function, n, x.get(), spline_c.get());
	}
	for(x_0 = show_a;x_0 - show_b < std::numeric_limits<double>::epsilon();x_0 += delta_x){
		if(mode != 3){
			y_0 = function(x_0);
			if(y_0 < min_y){
				min_y = y_0;
			}
			if(y_0 > max_y){
				max_y = y_0;
			}
		} else {
			if(n <= 50){
				y_0 = calculate_newton_discrepancy(x_0, function, n, x.get(), newton_c.get());
				if(y_0 < min_y){
					min_y = y_0;
				}
				if(y_0 > max_y){
					max_y = y_0;
				}
			}
			y_0 = calculate_spline_discrepancy(x_0, function, n, x.get(), spline_c.get());
			if(y_0 < min_y){
				min_y = y_0;
			}
			if(y_0 > max_y){
				max_y = y_0;
			}
		}
		if((mode == 0 || mode == 2) && n <= 50){
			y_0 = calculate_newton_approximation(x_0, n, x.get(), newton_c.get());
			if(y_0 < min_y){
				min_y = y_0;
			}
			if(y_0 > max_y){
				max_y = y_0;
			}
		}
		if(mode == 1 || mode == 2){
			y_0 = calculate_spline_approximation(x_0, n, x.get(), spline_c.get());
			if(y_0 < min_y){
				min_y = y_0;
			}
			if(y_0 > max_y){
				max_y = y_0;
			}
		}
	}
	
	max_abs_F = std::fmax(std::fabs(min_y), std::fabs(max_y));
	char str[24];
	std::snprintf(str, 24, "max{|F|} = %9.2e", max_abs_F);
	max_abs_F_changed(QString(str));
	
	double delta_y;
	if(std::fabs(max_y - min_y) <= std::numeric_limits<double>::epsilon() * (0.5 * std::fabs(max_y + min_y))){
		delta_y = 1e-14;
	} else {
		delta_y = 0.05 * (max_y - min_y);
	}
	min_y -= delta_y;
	max_y += delta_y;
}

void Window::init_status_bar(){
	function_description_changed(QString(get_function_description(k)));
	char str[24];
	std::snprintf(str, 16, "n = %d", n);
	number_of_points_changed(QString(str));
	std::snprintf(str, 16, "s = %d", s);
	scale_changed(QString(str));
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	std::snprintf(str, 24, "max{|F|} = %9.2e", max_abs_F);
	max_abs_F_changed(QString(str));
}

// change current function for drawing
void Window::change_function(){
	k++;
	if(k > 6){
		k = 0;
	}
	function = get_function(k);
	derivative = get_derivative(k);
	calculate_points();
	if(n <= 50){
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	}
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_max_abs_f();
	calculate_min_max();
	function_description_changed(QString(get_function_description(k)));
	p = 0;
	char str[16];
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	update();
}

void Window::change_graph(){
	mode++;
	if(mode >= 4){
		mode = 0;
	}
	calculate_min_max();
	update();
}

void Window::increase_scale(){
	s++;
	double center = (show_a + show_b) * 0.5;
	show_a = (show_a + center) * 0.5;
	show_b = (show_b + center) * 0.5;
	char str[16];
	std::snprintf(str, 16, "s = %d", s);
	scale_changed(QString(str));
	calculate_min_max();
	update();
}

void Window::decrease_scale(){
	s--;
	double center = (show_a + show_b) * 0.5;
	show_a = 2 * show_a - center;
	show_b = 2 * show_b - center;
	char str[16];
	std::snprintf(str, 16, "s = %d", s);
	scale_changed(QString(str));
	calculate_min_max();
	update();
}

void Window::increase_points(){
	n *= 2;
	if(n <= 0){
		n = 1;
	}
	x.reset(new double[n]);
	f.reset(new double[n]);
	d.reset(new double[n]);
	calculate_points();
	if(n <= 50){
		newton_c.reset(new double[2 * n]);
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	} else {
		newton_c.reset(nullptr);
	}
	spline_c.reset(new double[4 * n]);
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_min_max();
	char str[16];
	std::snprintf(str, 16, "n = %d", n);
	number_of_points_changed(QString(str));
	p = 0;
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	update();
}

void Window::decrease_points(){
	n /= 2;
	if(n <= 0){
		n = 1;
	}
	x.reset(new double[n]);
	f.reset(new double[n]);
	d.reset(new double[n]);
	calculate_points();
	if(n <= 50){
		newton_c.reset(new double[2 * n]);
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	} else {
		newton_c.reset(nullptr);
	}
	spline_c.reset(new double[4 * n]);
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_min_max();
	char str[16];
	std::snprintf(str, 16, "n = %d", n);
	number_of_points_changed(QString(str));
	p = 0;
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	update();
}

void Window::add_distortion(){
	p++;
	f[n / 2] += 0.1 * max_abs_f;
	char str[16];
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	if(n <= 50){
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	}
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_min_max();
	update();
}

void Window::subtract_distortion(){
	p--;
	f[n / 2] -= 0.1 * max_abs_f;
	char str[16];
	std::snprintf(str, 16, "p = %d", p);
	distortion_changed(QString(str));
	if(n <= 50){
		make_Lagrange_polynomial(n, x.get(), f.get(), d.get(), newton_c.get());
	}
	make_spline(n, x.get(), f.get(), spline_c.get());
	calculate_min_max();
	update();
}

QPointF Window::l2g(double x_0, double y_0){
	double x_gl = (x_0 - show_a) / (show_b - show_a) * width();
	double y_gl = (max_y - y_0) / (max_y - min_y) * height();
	return QPointF(x_gl, y_gl);
}

// draw approximated line for graph
template<typename Function>
void Window::drawGraph(QPainter& painter, Function func){
	double delta_x = (show_b - show_a) / width();
	double x1, x2, y1, y2;
	
	x1 = show_a;
	y1 = func(x1);
	for (x2 = x1 + delta_x; x2 - show_b < std::numeric_limits<double>::epsilon(); x2 += delta_x) {
		y2 = func(x2);
		// local coords are converted to draw coords
		painter.drawLine(l2g(x1, y1), l2g(x2, y2));

		x1 = x2;
		y1 = y2;
	}
	x2 = show_b;
	y2 = func(x2);
	painter.drawLine(l2g(x1, y1), l2g(x2, y2));
}

// render graph
void Window::paintEvent(QPaintEvent* event){
	(void)event;
	QPainter painter(this);
	QPen pen_black(Qt::black, 2);
	QPen pen_red(QColor("palevioletred"), 4);
	QPen pen_blue(QColor("royalblue"), 4);
	QPen pen_green(QColor("green"), 4);
	char str1[16];
	char str2[16];
	
	std::printf("k = %d %-20s n = %d p = %d s = %d [a, b] = [%10.3e, %10.3e] max{|F|} = %10.3e\n", k, get_function_description(k), n, p, s, show_a, show_b, max_abs_F);

	if(mode != 3){
		painter.setPen(pen_red);
		drawGraph(painter, function);
	} else {
		if(n <= 50){
			painter.setPen(pen_blue);
			drawGraph(painter, std::bind(calculate_newton_discrepancy, std::placeholders::_1, function, n, x.get(), newton_c.get()));
		}
		
		painter.setPen(pen_green);
		drawGraph(painter, std::bind(calculate_spline_discrepancy, std::placeholders::_1, function, n, x.get(), spline_c.get()));
	}
	if((mode == 0 || mode == 2) && n <= 50){
		painter.setPen(pen_blue);
		drawGraph(painter, std::bind(calculate_newton_approximation, std::placeholders::_1, n, x.get(), newton_c.get()));
	}
	if(mode == 1 || mode == 2){
		painter.setPen(pen_green);
		drawGraph(painter, std::bind(calculate_spline_approximation, std::placeholders::_1, n, x.get(), spline_c.get()));
	}

	// draw axis
	painter.setPen(pen_black);
	painter.drawLine(l2g(show_a, 0), l2g(show_b, 0));
	painter.drawLine(l2g(0, min_y), l2g(0, max_y));
	
	double show_a_label = show_a + 0.01 * (show_b - show_a);
	double show_b_label = show_b - 0.01 * (show_b - show_a);
	double min_y_label = min_y + 0.01 * (max_y - min_y);
	double max_y_label = max_y - 0.01 * (max_y - min_y);
	
	painter.drawLine(l2g(show_a_label, -0.02 * (max_y - min_y)), l2g(show_a_label, 0.02 * (max_y - min_y)));
	painter.drawLine(l2g(show_b_label, -0.02 * (max_y - min_y)), l2g(show_b_label, 0.02 * (max_y - min_y)));
	
	painter.drawLine(l2g(-0.02 * (show_b - show_a), min_y_label), l2g(0.02 * (show_b - show_a), min_y_label));
	painter.drawLine(l2g(-0.02 * (show_b - show_a), max_y_label), l2g(0.02 * (show_b - show_a), max_y_label));
	
	if((show_b - show_a) < 1e-2 || (show_b - show_a) > 1e2){
		std::snprintf(str1, 16, "%9.2e", show_a_label);
		std::snprintf(str2, 16, "%9.2e", show_b_label);
	} else {
		std::snprintf(str1, 16, "%.2f", show_a_label);
		std::snprintf(str2, 16, "%.2f", show_b_label);
	}
	if(l2g(show_a_label, 0.03 * (max_y - min_y)).y() + (painter.fontMetrics().height() / 2) >= 0 && l2g(show_a_label, 0.03 * (max_y - min_y)).y() <= height()){
		painter.drawText(l2g(show_a_label, 0.03 * (max_y - min_y)), str1);
		painter.drawText(l2g(show_b_label, 0.03 * (max_y - min_y)) - QPointF(painter.fontMetrics().boundingRect(str2).width(), 0), str2);
	}
	
	if((max_y - min_y) < 1e-2 || (max_y - min_y) > 1e2){
		std::snprintf(str1, 16, "%9.2e", min_y_label);
		std::snprintf(str2, 16, "%9.2e", max_y_label);
	} else {
		std::snprintf(str1, 16, "%.2f", min_y_label);
		std::snprintf(str2, 16, "%.2f", max_y_label);
	}
	if(l2g(0.03 * (show_b - show_a), min_y_label).x() - painter.fontMetrics().boundingRect(str2).width() >= 0 && l2g(0.03 * (show_b - show_a), min_y_label).x() <= width()){
		painter.drawText(l2g(0.03 * (show_b - show_a), min_y_label), str1);
		painter.drawText(l2g(0.03 * (show_b - show_a), max_y_label) + QPointF(0, painter.fontMetrics().height() / 2), str2);
	}
}
