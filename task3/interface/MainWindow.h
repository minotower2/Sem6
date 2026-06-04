#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QHBoxLayout>

#include "utils/CommonDefs.h"
#include "gui/GraphWidget.h"
#include "gui/InfoDisplay.h"
#include "model/Model.h"
#include "resources/Storage.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT
    public:
        explicit MainWindow(double, double, double, double, int, int, int, int, int, QWidget *parent = nullptr);
        ~MainWindow();

    private:
		QWidget *central;
		GraphWidget *GUI;
		Model *model;
		Storage *storage;
		InfoDisplay *infoDisplay;
		Data *data;
};

#endif // MAINWINDOW_H
