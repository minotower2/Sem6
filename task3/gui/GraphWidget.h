#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QWidget>
#include <QPainter>
#include <QColor>
#include <QResizeEvent>
#include <QPolygonF>

#include <iostream>
//#include <cmath>

#include "gui/3DRender.h"
#include "gui/Data.h"

#include "utils/CommonDefs.h"
#include "utils/Geometry.h"
#include "resources/Storage.h"

class GraphWidget : public QWidget 
{
	Q_OBJECT
private:
	using Polygon = Storage::Polygon;

	Storage *storage;
	Data *data;

	QColor bgColor;

	int w, h;
	double dist;

	R3Geometry::Point center;
	Camera camera;

	R3Geometry::Point lightPos;
	PointLight light;

	bool rotating = false;
	QPoint lastMousePos;

	Polygon pointToPixels(const R3Geometry::Point &) const;
	Polygon segToPixels(const R3Geometry::Point &, const R3Geometry::Point &) const;
	Polygon triangleToPixels(const R3Geometry::Point &, const R3Geometry::Point &, const R3Geometry::Point &) const;

public:
	GraphWidget(QWidget *, double, double, double, double, Storage *, Data *);

	void scaleUp();
	void scaleDown();

	void rotateLeft();
	void rotateRight();
	
	QSize minimumSizeHint() const;
	QSize sizeHint() const;

protected:
	void paintEvent(QPaintEvent *);
	void resizeEvent(QResizeEvent *);

private:
	void drawCoordSystem(QPainter *);
	
	void drawMesh(QPainter *, const std::vector<Storage::Triangle> &, const std::vector<R3Geometry::Point> &);

	void mousePressEvent(QMouseEvent *);
	void mouseReleaseEvent(QMouseEvent *);
	void mouseMoveEvent(QMouseEvent *);

	void wheelEvent(QWheelEvent *);

	void initializeMap();

	QColor converstion(const PointLight::Color &) const;
private slots:
	void updated(void) 
	{
		update();
	}
};

#endif // GRAPHWIDGET_H
