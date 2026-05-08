#include <memory>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "window.h"

int main (int argc, char** argv){
	double a, b;
	int n, k;
	if(!(argc == 5
		&& sscanf(argv[1], "%lf", &a) == 1
		&& sscanf(argv[2], "%lf", &b) == 1 && a < b
		&& sscanf(argv[3], "%d", &n) == 1 && n > 0
		&& sscanf(argv[4], "%d", &k) == 1 && k >= 0 && k <= 6
	)){
		std::printf("Usage: %s (double a) (double b > a) (int n > 0) (int k in [0; 6])\n", argv[0]);
		return 1;
	}
	
	QApplication app(argc, argv);

	std::unique_ptr<QMainWindow> window = std::make_unique<QMainWindow>();
	QMenuBar* tool_bar = new QMenuBar(window.get());
	Window* graph_area = new Window(window.get(), a, b, n, k);
	QStatusBar* status_bar = new QStatusBar(window.get());
	QLabel* function_description = new QLabel;
	QLabel* number_of_points = new QLabel;
	QLabel* distortion = new QLabel;
	QLabel* scale = new QLabel;
	QLabel* max_abs_F = new QLabel;
	QAction *action;

	action = tool_bar->addAction("Change function", graph_area, SLOT(change_function()));
	action->setShortcut(QString("0"));
	action = tool_bar->addAction("Change graph", graph_area, SLOT(change_graph()));
	action->setShortcut(QString("1"));
	action = tool_bar->addAction("Double scale", graph_area, SLOT(increase_scale()));
	action->setShortcut(QString("2"));
	action = tool_bar->addAction("Half scale", graph_area, SLOT(decrease_scale()));
	action->setShortcut(QString("3"));
	action = tool_bar->addAction("Double approximation points", graph_area, SLOT(increase_points()));
	action->setShortcut(QString("4"));
	action = tool_bar->addAction("Half approximation points", graph_area, SLOT(decrease_points()));
	action->setShortcut(QString("5"));
	action = tool_bar->addAction("Add distortion", graph_area, SLOT(add_distortion()));
	action->setShortcut(QString("6"));
	action = tool_bar->addAction("Subtract distortion", graph_area, SLOT(subtract_distortion()));
	action->setShortcut(QString("7"));
	
	action = tool_bar->addAction("Exit", window.get(), SLOT(close()));
	action->setShortcut(QString("esc"));
	
	status_bar->insertPermanentWidget(0, function_description, 3);
	QObject::connect(graph_area, SIGNAL(function_description_changed(const QString&)), function_description, SLOT(setText(const QString&)));
	status_bar->insertPermanentWidget(1, number_of_points, 1);
	QObject::connect(graph_area, SIGNAL(number_of_points_changed(const QString&)), number_of_points, SLOT(setText(const QString&)));
	status_bar->insertPermanentWidget(2, scale, 1);
	QObject::connect(graph_area, SIGNAL(scale_changed(const QString&)), scale, SLOT(setText(const QString&)));
	status_bar->insertPermanentWidget(3, distortion, 1);
	QObject::connect(graph_area, SIGNAL(distortion_changed(const QString&)), distortion, SLOT(setText(const QString&)));
	status_bar->insertPermanentWidget(4, max_abs_F, 0);
	QObject::connect(graph_area, SIGNAL(max_abs_F_changed(const QString&)), max_abs_F, SLOT(setText(const QString&)));

	//tool_bar->setMaximumHeight(30);
	status_bar->setStyleSheet("background-color: rgb(224, 224, 224);");
	graph_area->init_status_bar();

	window->setMenuBar(tool_bar);
	window->setStatusBar(status_bar);
	window->setCentralWidget(graph_area);
	window->setWindowTitle("Graph");

	window->show();
	app.exec();
	return 0;
}
