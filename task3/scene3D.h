#ifndef SCENE3D_H // Защита от двойного включения заголовочного файла
#define SCENE3D_H //

#include <QGLWidget>     // Базовый класс для виджетов OpenGL в Qt
#include <QMouseEvent>   // Для обработки событий мыши
#include <QKeyEvent>     // Для обработки событий клавиатуры
#include <QWheelEvent>   // Для обработки событий колеса мыши
#include <QPoint>        // Для представления точки (координат)
#include <vector>        // Для std::vector в solveTridiagonalSystem

#define MIN_SPLINE_NODES 3 // Минимальное количество узлов для сплайна с экстраполяцией по 3 точкам

class Scene3D : public QGLWidget // Объявление класса Scene3D, наследуемого от QGLWidget
{
    Q_OBJECT // Макрос Qt, необходимый для использования сигналов и слотов

private:
    // Параметры сцены
    int nx_param;       // Желаемое количество узлов по оси X для интерполяции
    int ny_param;       // Желаемое количество узлов по оси Y для интерполяции
    int id;             // Идентификатор текущей отображаемой функции (0-7)
    int scale;          // Масштаб, управляемый пользователем для отрисовки (для области a,b,c,d)
    const char *f_name; // Имя текущей функции (например, "f(x,y) = 1")
    double eps_param;   // Параметр точности (эпсилон), считываемый из командной строки
    double a;           // Минимальное значение X для области построения функции
    double b;           // Максимальное значение X для области построения функции
    double c;           // Минимальное значение Y для области построения функции
    double d;           // Максимальное значение Y для области построения функции
    int p = 0;          // Параметр для внесения искусственной ошибки (выброса) в центре графика
    double fmax_orig;   // Максимальное абсолютное значение исходной функции на заданной области
    int gr = 0;         // Тип отображаемого графика (0: исходная функция, 1: аппроксимация, 2: ошибка)
    int repeat;         // Угол поворота для отображения (связан с дискретным поворотом)
    int s = 0;          // Параметр для отображения масштабирования nx/ny (влияет на текст)
    double max_metod;   // Максимальное абсолютное значение аппроксимированной функции

    // Указатель на функцию
    double (*f)(double, double);     // Указатель на текущую функцию f(x,y)

    // Параметры OpenGL
    GLfloat xRot; // Угол поворота вокруг оси X
    GLfloat yRot; // Угол поворота вокруг оси Y
    GLfloat zRot; // Угол поворота вокруг оси Z
    GLfloat zTra; // Смещение вдоль оси Y (в OpenGL координатах это часто ось Z сцены, но здесь используется для Y-смещения камеры)
    GLfloat nSca; // Фактор масштабирования OpenGL

    QPoint ptrMousePosition; // Позиция указателя мыши для отслеживания перемещений

    // Данные для интерполяции сплайнами (бикубической Эрмитовой)
    int current_nx;     // Фактическое количество узлов по X в сетке сплайна
    int current_ny;     // Фактическое количество узлов по Y в сетке сплайна
    double* spline_x_nodes; // Массив X-координат узлов сетки
    double* spline_y_nodes; // Массив Y-координат узлов сетки
    double** f_val;     // Значения функции f(xi, yj) в узлах сетки
    double** fx_val;    // Значения производной df/dx (xi, yj) в узлах сетки (из 1D сплайнов)
    double** fy_val;    // Значения производной df/dy (xi, yj) в узлах сетки (из 1D сплайнов)
    double** fxy_val;   // Значения смешанной производной d2f/dxdy (xi, yj) в узлах сетки (из 1D сплайнов)

    // Вспомогательные методы для сплайнов
    void cleanupSplineData(); // Освобождение памяти, выделенной под данные сплайна
    void allocateSplineData(int rows, int cols); // Выделение памяти для данных сплайна
    // Вычисляет 1D кубический сплайн и возвращает ПЕРВЫЕ производные в узлах.
    // Граничные условия (вторые производные M0, Mn-1) вычисляются внутри экстраполяцией.
    void compute_1D_Spline_derivatives(const double* input_vals, const double* coords, int N,
                                       double* output_first_derivatives);
    // Решает СЛАУ Ax=B для трехдиагональной матрицы методом прогонки.
    // a - поддиагональ (a[0] не используется), d - главная диагональ, c - наддиагональ, b - правая часть, x_sol - решение, n_sys - размер системы.
    bool solveTridiagonalSystem(const double* a, double* d_mod, const double* c, double* b_mod,
                                double* x_sol, int n_sys);
    // Экстраполирует вторую производную на границе (M0 или Mn-1) используя 3 точки.
    double extrapolate_M_boundary(const double* coords, const double* values, int N_total, bool at_start_node);

    void computeSplineData();    // Вычисление всех необходимых данных для интерполяции (f_val, fx_val, fy_val, fxy_val)
    double evaluateSpline(double x, double y); // Вычисление значения аппроксимирующей функции в точке (x,y)
    double errorSpline(double x, double y); // Вычисление абсолютной ошибки между аппроксимацией и исходной функцией

    // Методы управления сценой
    void scale_plus();      // Увеличить масштаб отображения (OpenGL nSca)
    void scale_minus();     // Уменьшить масштаб отображения (OpenGL nSca)
    void rotate_up();       // Повернуть сцену вверх
    void rotate_down();     // Повернуть сцену вниз
    void rotate_left();     // Повернуть сцену влево
    void rotate_right();    // Повернуть сцену вправо
    void p_angle();         // Увеличить угол дискретного поворота (zRot)
    void m_angle();         // Уменьшить угол дискретного поворота (zRot)
    void defaultScene();    // Сбросить параметры вида сцены по умолчанию
    void change_func();     // Сменить отображаемую функцию (циклически)
    void change_graph();    // Сменить тип отображаемого графика (исходная, аппроксимация, ошибка)
    void increase_param_scale(); // Увеличить масштаб параметров (сузить область a,b,c,d)
    void decrease_param_scale(); // Уменьшить масштаб параметров (расширить область a,b,c,d)
    void increase_nx_ny_param(); // Увеличить количество узлов сетки (nx_param, ny_param)
    void decrease_nx_ny_param(); // Уменьшить количество узлов сетки (nx_param, ny_param)
    void increase_error_p();     // Увеличить параметр искусственной ошибки 'p'
    void decrease_error_p();     // Уменьшить параметр искусственной ошибки 'p'
    void plus_angle();           // Обертка для p_angle с обновлением GL
    void minus_angle();          // Обертка для m_angle с обновлением GL


    // Методы отрисовки
    void drawAxis();                // Отрисовка осей координат
    void drawOriginalFunction();    // Отрисовка исходной функции
    void drawSplineApproximation(); // Отрисовка аппроксимации сплайнами
    void drawErrorFunction();       // Отрисовка функции ошибки

protected:
    // Переопределенные методы QGLWidget
    void initializeGL() override; // Инициализация OpenGL состояния
    void resizeGL(int nWidth, int nHeight) override; // Обработка изменения размеров виджета
    void paintGL() override;      // Основной метод отрисовки сцены

    // Обработчики событий
    void mousePressEvent(QMouseEvent *pe) override;   // Нажатие кнопки мыши
    void mouseMoveEvent(QMouseEvent *pe) override;    // Перемещение мыши с зажатой кнопкой
    void mouseReleaseEvent(QMouseEvent *pe) override; // Отпускание кнопки мыши
    void wheelEvent(QWheelEvent *pe) override;        // Прокрутка колеса мыши
    void keyPressEvent(QKeyEvent *pe) override;       // Нажатие клавиши на клавиатуре

public:
    double max_original_func(); // Вычисляет максимальное абсолютное значение исходной функции
    ~Scene3D();                 // Деструктор класса
    Scene3D(char *argv[], QWidget *parent = nullptr); // Конструктор класса
    int parse_command_line(int argc, char **argv); // Разбор аргументов командной строки
};

#endif // SCENE3D_H
