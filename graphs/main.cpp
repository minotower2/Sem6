 
#include <memory> // Для std::unique_ptr
#include <QtWidgets/QApplication> // Класс приложения Qt
#include <QtWidgets/QMainWindow> // Главное окно
#include <QtWidgets/QVBoxLayout> // Компоновщик (хотя не используется напрямую в этом коде, т.к. Window установлен как CentralWidget)
#include <QtWidgets/QAction> // Действия (пункты меню, сочетания клавиш)
#include <QtWidgets/QMenuBar> // Строка меню
#include <QtWidgets/QMessageBox> // Диалоговые окна сообщений (пока не используются, но класс подключен)
#include <QtWidgets/QLabel> // Текстовые метки для статусбара
#include <QtWidgets/QStatusBar> // Строка состояния
#include <cstdio> // Для sscanf, printf
#include "window.h" // Заголовочный файл нашего главного виджета

int main (int argc, char** argv){ // Главная функция программы
	double a, b; // Переменные для границ отрезка
	int n, k; // Переменные для числа точек и номера функции

    // Проверяем количество аргументов командной строки и парсим их
    // argc - количество аргументов (имя программы + 4 параметра)
    // argv - массив строк аргументов
	if(!(argc == 5 // Должно быть ровно 5 аргументов
		&& sscanf(argv[1], "%lf", &a) == 1 // Читаем первый аргумент как double в 'a'
		&& sscanf(argv[2], "%lf", &b) == 1 && a < b // Читаем второй аргумент как double в 'b', проверяем, что b > a
		&& sscanf(argv[3], "%d", &n) == 1 && n > 0 // Читаем третий аргумент как int в 'n', проверяем, что n > 0
		&& sscanf(argv[4], "%d", &k) == 1 && k >= 0 && k <= 6 // Читаем четвертый аргумент как int в 'k', проверяем, что k в диапазоне [0, 6]
	)){
		// Если аргументы некорректны, выводим сообщение об использовании и завершаем программу
		std::printf("Надо так: %s (double a) (double b > a) (int n > 0) (int k в [0; 6])\n", argv[0]);
		return 1; // Код ошибки
	}

	QApplication app(argc, argv); // Создаем объект приложения Qt. Он управляет ресурсами приложения и циклом событий.

	std::unique_ptr<QMainWindow> window = std::make_unique<QMainWindow>(); // Создаем главное окно приложения (с помощью умного указателя)
	QMenuBar* tool_bar = new QMenuBar(window.get()); // Создаем строку меню, родителем указываем главное окно
	Window* graph_area = new Window(window.get(), a, b, n, k); // Создаем наш пользовательский виджет Window, передавая ему параметры и главное окно как родителя
	QStatusBar* status_bar = new QStatusBar(window.get()); // Создаем строку состояния, родителем указываем главное окно

    // Создаем QLabel виджеты, которые будут отображаться в строке состояния
	QLabel* function_description = new QLabel; // Для описания функции
	QLabel* number_of_points = new QLabel; // Для числа точек n
	QLabel* distortion = new QLabel; // Для параметра искажения p
	QLabel* scale = new QLabel; // Для параметра масштаба s
	QLabel* max_abs_F = new QLabel; // Для максимальной абсолютной невязки/значения

	QAction *action; // Указатель на действие (используется временно для создания действий)

    // Создаем действия (Action) для меню и сочетаний клавиш, связываем их с методами (слотами) нашего graph_area
	action = tool_bar->addAction("Change function", graph_area, SLOT(change_function())); // Действие "Change function", вызывает слот change_function(), родитель - tool_bar, получатель слота - graph_area
	action->setShortcut(QString("0")); // Устанавливаем сочетание клавиш '0'

	action = tool_bar->addAction("Change graph", graph_area, SLOT(change_graph())); // Действие "Change graph"
	action->setShortcut(QString("1")); // Клавиша '1'

	action = tool_bar->addAction("Double scale", graph_area, SLOT(increase_scale())); // Действие "Double scale"
	action->setShortcut(QString("2")); // Клавиша '2'

	action = tool_bar->addAction("Half scale", graph_area, SLOT(decrease_scale())); // Действие "Half scale"
	action->setShortcut(QString("3")); // Клавиша '3'

	action = tool_bar->addAction("Double approximation points", graph_area, SLOT(increase_points())); // Действие "Double points"
	action->setShortcut(QString("4")); // Клавиша '4'

	action = tool_bar->addAction("Half approximation points", graph_area, SLOT(decrease_points())); // Действие "Half points"
	action->setShortcut(QString("5")); // Клавиша '5'

	action = tool_bar->addAction("Add distortion", graph_area, SLOT(add_distortion())); // Действие "Add distortion"
	action->setShortcut(QString("6")); // Клавиша '6'

	action = tool_bar->addAction("Subtract distortion", graph_area, SLOT(subtract_distortion())); // Действие "Subtract distortion"
	action->setShortcut(QString("7")); // Клавиша '7'

	action = tool_bar->addAction("Exit", window.get(), SLOT(close())); // Действие "Exit"
	action->setShortcut(QString("esc")); // Клавиша Escape. Закрывает главное окно.

    // Вставляем QLabel виджеты в строку состояния
    // insertPermanentWidget(индекс, виджет, растягивающий фактор) - вставляет виджет слева направо как постоянный элемент
	status_bar->insertPermanentWidget(0, function_description, 3); // Описание функции (растягивается в 3 раза больше, чем следующий)
    // Соединяем сигнал от graph_area function_description_changed() со слотом setText() QLabel виджета.
	QObject::connect(graph_area, SIGNAL(function_description_changed(const QString&)), function_description, SLOT(setText(const QString&)));

	status_bar->insertPermanentWidget(1, number_of_points, 1); // n
	QObject::connect(graph_area, SIGNAL(number_of_points_changed(const QString&)), number_of_points, SLOT(setText(const QString&)));

	status_bar->insertPermanentWidget(2, scale, 1); // s
	QObject::connect(graph_area, SIGNAL(scale_changed(const QString&)), scale, SLOT(setText(const QString&)));

	status_bar->insertPermanentWidget(3, distortion, 1); // p
	QObject::connect(graph_area, SIGNAL(distortion_changed(const QString&)), distortion, SLOT(setText(const QString&)));

	status_bar->insertPermanentWidget(4, max_abs_F, 0); // max{|F|} (не растягивается, занимает минимум места)
	QObject::connect(graph_area, SIGNAL(max_abs_F_changed(const QString&)), max_abs_F, SLOT(setText(const QString&)));

	//tool_bar->setMaximumHeight(30); // Устанавливаем максимальную высоту для панели инструментов (закомментировано)
	status_bar->setStyleSheet("background-color: rgb(224, 224, 224);"); // Задаем стиль (цвет фона) для строки состояния

	graph_area->init_status_bar(); // Инициализируем текст в строке состояния (вызываем метод Window)

	window->setMenuBar(tool_bar); // Устанавливаем строку меню в главное окно
	window->setStatusBar(status_bar); // Устанавливаем строку состояния в главное окно
	window->setCentralWidget(graph_area); // Устанавливаем наш виджет graph_area как центральный виджет окна
	window->setWindowTitle("Graph"); // Устанавливаем заголовок окна

	window->show(); // Отображаем главное окно со всеми его содержимым

	app.exec(); // Запускаем главный цикл обработки событий Qt. Программа будет работать до закрытия окна.
	return 0; // Возвращаем 0 при успешном завершении
}
