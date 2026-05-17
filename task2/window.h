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
	double max_abs_f = 0; // Макс. абс. значение ИСХОДНОЙ функции на [a,b]
	double max_abs_F = 0; // Макс. абс. значение ОТОБРАЖАЕМЫХ функций/невязок на [show_a, show_b]

	std::unique_ptr<double[]> x = nullptr; // Узлы для Эрмита
	std::unique_ptr<double[]> f = nullptr; // Значения f(x_i) для Эрмита
	std::unique_ptr<double[]> d = nullptr; // Значения d(x_i) для Эрмита

    std::unique_ptr<double[]> chebyshev_c = nullptr; // Коэффициенты для Чебышева
    std::unique_ptr<double[]> hermite_c = nullptr;    // Коэффициенты для сплайна Бесселя

	double (*function)(double); // Указатель на исходную функцию
	double (*derivative)(double); // Указатель на производную

	int n = 10; // Количество узлов (для Эрмита и как параметр для Чебышева)
	int k = 0;  // Индекс функции
	int s = 0;  // Масштаб
	int p = 0;  // Искажение (влияет на f.get() для Эрмита)
    
    // Новые режимы:
    // 0: Исходная функция
    // 1: Исходная + Чебышев
    // 2: Исходная + Эрмит
    // 3: Исходная + Чебышев + Эрмит
    // 4: Невязка Чебышева
    // 5: Невязка Эрмита
    // 6: Обе невязки
	int mode = 0; 
    int total_modes = 7; // Количество режимов от 0 до 6

public:
	Window(QWidget *parent, double a, double b, int n, int k);

	QSize minimumSizeHint() const {
		return QSize(100, 100);
	}
	QSize sizeHint() const {
		return QSize(1000, 1000);
	}

	void calculate_points(); // Для узлов Бесселя
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
