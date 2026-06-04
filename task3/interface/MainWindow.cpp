#include "MainWindow.h"
#include <locale.h>

MainWindow::MainWindow(double ax, double bx, double ay, double by, int nx, int ny, int mx, int my, int k, QWidget *parent): QMainWindow(parent)
{
    setlocale(LC_ALL, "C");
	
	storage = new Storage(this);
	data = new Data(this);

	model = new Model(ax, bx, ay, by, nx, ny, mx, my, k, storage, data, this);

	central = new QWidget(this);
	GUI = new GraphWidget (this, ax, bx, ay, by, storage, data);
	infoDisplay = new InfoDisplay(this, data);

	infoDisplay->setFixedWidth(250);

	QHBoxLayout *layout = new QHBoxLayout(central);

	layout->setSpacing(0);
	layout->setContentsMargins(0, 0, 0, 0);

	layout->addWidget(GUI);
	layout->addWidget(infoDisplay);

	setCentralWidget (central);

	QMenu *fileMenu = menuBar()->addMenu("&File");
	QMenu *funcMenu = menuBar()->addMenu("&Function");

	QAction *changeFuncAction = funcMenu->addAction("Change function", model, &Model::changeFunc);
	changeFuncAction->setShortcut(QKeySequence("Ctrl+0"));

	QAction *changeTypeAction = funcMenu->addAction("Change type", model, &Model::changeMode);
	changeTypeAction->setShortcut(QKeySequence("Ctrl+1"));

	QAction *scaleUpAction = funcMenu->addAction("Scale up", GUI, &GraphWidget::scaleUp);
	scaleUpAction->setShortcut(QKeySequence("Ctrl+2"));

	QAction *scaleDownAction = funcMenu->addAction("Scale down", GUI, &GraphWidget::scaleDown);
	scaleDownAction->setShortcut(QKeySequence("Ctrl+3"));

	QAction *numberUpAction = funcMenu->addAction("Number of points up", model, &Model::numberUp);
	numberUpAction->setShortcut(QKeySequence("Ctrl+4"));

	QAction *numberDownAction = funcMenu->addAction("Number of points down", model, &Model::numberDown);
	numberDownAction->setShortcut(QKeySequence("Ctrl+5"));
	
	QAction *triangulationUpAction = funcMenu->addAction("Number of points of triangulation up", model, &Model::upTriangulation);
	triangulationUpAction->setShortcut(QKeySequence("Ctrl+U"));

	QAction *triangulationDownAction = funcMenu->addAction("Number of points of triangulation down", model, &Model::downTriangulation);
	triangulationDownAction->setShortcut(QKeySequence("Ctrl+D"));

	QAction *errorUpAction = funcMenu->addAction("Error up", model, &Model::errorUp);
	errorUpAction->setShortcut(QKeySequence("Ctrl+6"));

	QAction *errorDownAction = funcMenu->addAction("Error down", model, &Model::errorDown);
	errorDownAction->setShortcut(QKeySequence("Ctrl+7"));

	QAction *rotateLeftAction = funcMenu->addAction("Rotate left on Pi/12", GUI, &GraphWidget::rotateLeft);
	rotateLeftAction->setShortcut(QKeySequence("Ctrl+8"));

	QAction *rotateRightAction = funcMenu->addAction("Rotate right on Pi/12", GUI, &GraphWidget::rotateRight);
	rotateRightAction->setShortcut(QKeySequence("Ctrl+9"));

	QAction *exitAction = fileMenu->addAction("E&xit", this, &QWidget::close);
	exitAction->setShortcut(QKeySequence::Quit);

	setWindowTitle ("Graph");
}

MainWindow::~MainWindow()
{
}
