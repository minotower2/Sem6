#ifndef WINDOW_H
#define WINDOW_H

#include <memory>
#include <QtWidgets/QtWidgets>
#include <QLabel>

class Window : public QWidget {
	Q_OBJECT

private:
	double a = -1.;
	double show_a = -1.;
	double b = 1.;
	double show_b = 1.;
	double min_y = 0;
	double max_y = 0;
	double max_abs_f = 0;
	double max_abs_F = 0;
	std::unique_ptr<double[]> x = nullptr;
	std::unique_ptr<double[]> f = nullptr;
	std::unique_ptr<double[]> d = nullptr;
	std::unique_ptr<double[]> newton_c = nullptr;
	std::unique_ptr<double[]> spline_c = nullptr;
	double (*function)(double);
	double (*derivative)(double);
	int n = 10;
	int k = 0;
	int s = 0;
	int p = 0;
	int mode = 0;

public:
	Window(QWidget *parent, double a, double b, int n, int k);

	QSize minimumSizeHint() const {
		return QSize(100, 100);
	}
	QSize sizeHint() const {
		return QSize(1000, 1000);
	}

	void calculate_points();
	void calculate_max_abs_f();
	void calculate_min_max();
	void init_status_bar();
public slots:
	void change_function();
	void change_graph();
	void increase_scale();
	void decrease_scale();
	void increase_points();
	void decrease_points();
	void add_distortion();
	void subtract_distortion();
signals:
	void function_description_changed(const QString&);
	void number_of_points_changed(const QString&);
	void distortion_changed(const QString&);
	void scale_changed(const QString&);
	void max_abs_F_changed(const QString&);

protected:
	QPointF l2g(double x_0, double y_0);
	template<typename Function>
	void drawGraph(QPainter& painter, Function func);
	void paintEvent(QPaintEvent *event);
};

#endif
