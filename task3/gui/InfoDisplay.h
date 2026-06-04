#ifndef INFODISPLAY_H
#define INFODISPLAY_H

#include <QWidget>
#include <QPainter>
#include <QColor>
#include <QPaintEvent>
#include <QVBoxLayout>
#include <QLabel>

#include <iostream>
#include "gui/Data.h"

class InfoDisplay : public QWidget 
{
	Q_OBJECT
private:
	Data *data;
	QColor bgColor;

	QVBoxLayout *table;
	QLabel *label;

public:
	InfoDisplay(QWidget *, Data *);
	
protected:
	void paintEvent(QPaintEvent *);

private:
	void drawInformation(void);
	void drawConsole(void);

public slots:
	void dataChanged() {
		update();
	}
};
#endif
