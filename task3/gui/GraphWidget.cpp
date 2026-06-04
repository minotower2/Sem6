#include "GraphWidget.h"

GraphWidget::GraphWidget (QWidget *parent, double ax, double bx, double ay, double by, Storage *storage, Data *data) :
	QWidget(parent),
	storage(storage), data(data),
	bgColor(Qt::white),
	dist(2*sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay) + 4)),
	center((ax+bx)/2, (ay+by)/2, 0),
	camera(center, Quaternion(1, 0, 0, 0), 1, dist, 2.0),
	lightPos({0, 0, 3}),
	light(lightPos)
{
	setAttribute(Qt::WA_StaticContents); // for optimizing painting events
	connect(storage, &Storage::updated, this, &GraphWidget::updated);
	w = width();
	h = height();
	camera.setAspect(w * (1.0)/h);
	camera.orbitPitch(-3*M_PI/4);
	camera.orbitYaw(3*M_PI/4);
}

void GraphWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);

	painter.setRenderHint(QPainter::Antialiasing, true);
	initializeMap();
	drawCoordSystem(&painter);

	auto mesh = storage->renderrAccess();
	
	drawMesh(&painter, mesh->triangles, mesh->points);
	painter.drawText(20, 20, QString::fromStdString(camera.orien()));	
	event->accept();
}


void GraphWidget::drawCoordSystem(QPainter *painter) 
{
	painter->setBrush(QBrush(bgColor));
	painter->drawRect(0, 0, w, h);
	painter->setBrush(Qt::NoBrush);

	double L = 10*dist;

	QPen penBlack(Qt::black);
	penBlack.setWidth(1);
	painter->setPen(penBlack);
	
	auto OX = segToPixels({-L, 0, 0}, {L, 0, 0});
	if (OX.visible)
		painter->drawLine(OX.vertex[0], OX.vertex[1]);

	auto OY = segToPixels({0, -L, 0}, {0, L, 0});
	if (OY.visible)
		painter->drawLine(OY.vertex[0], OY.vertex[1]);

	auto OZ = segToPixels({0, 0, -L}, {0, 0, L});
	if (OZ.visible)
		painter->drawLine(OZ.vertex[0], OZ.vertex[1]);

	auto zero = pointToPixels({0, 0, 0});
	if (zero.visible)
		painter->drawText(zero.vertex[0].x()+10, zero.vertex[0].y()-10, QString("0"));

	if (OX.visible)
		painter->drawText(OX.vertex[1]+QPointF((w/2 - OX.vertex[1].x() < 0 ? -1 : 1)*15, -10), QString("X"));

	if (OY.visible)
		painter->drawText(OY.vertex[0]+QPointF((w/2 - OY.vertex[0].x() < 0 ? -1 : 1)*15, -10), QString("Y"));

	if (OZ.visible)
		painter->drawText(OZ.vertex[1]+QPointF(10, (h/2 - OZ.vertex[1].y() < 0 ? -1 : 1)*15), QString("Z"));
	/*

	QPointF center = mapToPixels(Approx::Point(0, 0));
	
	QPointF dx(8., 0.);
	QPointF dy(0., 8.);

	QPointF xTick(100., 0);
	QPointF yTick(0, 100.);

	int i = 0;
	while ((center+(++i)*xTick).x() < w) {
		painter->drawLine((center+i*xTick)-dy, (center+i*xTick)+dy);
		painter->drawText((center+i*xTick).x(), (center+i*xTick).y() - 20, QString::fromStdString(std::to_string(mapFromPixels(center+i*xTick).X())));
	}

	i = 0;
	while ((center+(--i)*xTick).x() > 0) {
		painter->drawLine((center+i*xTick)-dy, (center+i*xTick)+dy);
		painter->drawText((center+i*xTick).x(), (center+i*xTick).y() - 20, QString::fromStdString(std::to_string(mapFromPixels(center+i*xTick).X())));
	}

	i = 0;
	while ((center+(++i)*yTick).y() < h) {
		painter->drawLine((center+i*yTick)-dx, (center+i*yTick)+dx);
		painter->drawText((center+i*yTick).x() + 20, (center+i*yTick).y(), QString::fromStdString(std::to_string(mapFromPixels(center+i*yTick).Y())));
	}

	i = 0;
	while ((center+(--i)*yTick).y() > 0) {
		painter->drawLine((center+i*yTick)-dx, (center+i*yTick)+dx);
		painter->drawText((center+i*yTick).x() + 20, (center+i*yTick).y(), QString::fromStdString(std::to_string(mapFromPixels(center+i*yTick).Y())));
	}
	*/
}

void GraphWidget::drawMesh(QPainter *painter, const std::vector<Storage::Triangle> &triangles, const std::vector<R3Geometry::Point> &points)
{
	QColor currColor;

	auto polygons = storage->poly();
	polygons.clear();

	for (const auto& tr : triangles) {
		polygons.push_back(triangleToPixels(points[tr.a], points[tr.b], points[tr.c]));
	}

	std::sort(polygons.begin(), polygons.end(), [](const Storage::Polygon &a, const Storage::Polygon &b) {return a.depth < b.depth;});

	for (const auto& poly : polygons) {
		if (!poly.visible) 
			continue;
		
		if (poly.color != currColor) {
			painter->setBrush(poly.color);
			painter->setPen(poly.color);
			currColor = poly.color;
		}
		painter->drawPolygon(poly.vertex, poly.n);
	}
}

void GraphWidget::resizeEvent(QResizeEvent *event) 
{
	if (width() == 0 || height() == 0)
		return;
		
	initializeMap();
	update();	
	event->accept();
}

void GraphWidget::initializeMap() 
{
	w = width();
	h = height();
	camera.setAspect(w *(1.)/h);
}

void GraphWidget::scaleUp()
{
	if (dist <= 25) {
		dist *= 2;
		camera.setDist(dist);
	}

	update();
}

void GraphWidget::scaleDown()
{
	if (dist >= 4) {
		dist /= 2;
		camera.setDist(dist);
	}

	update();
}

void GraphWidget::rotateLeft()
{
	camera.orbitYaw(M_PI/12);
	update();
}

void GraphWidget::rotateRight()
{
	camera.orbitYaw(-M_PI/12);
	update();
}


void GraphWidget::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) {
		rotating = true;
		lastMousePos = event->pos();
	}
}

void GraphWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton)
		rotating = false;
}

void GraphWidget::mouseMoveEvent(QMouseEvent *event)
{
	if (!rotating)
		return;

	QPoint delta = event->pos() - lastMousePos;
	lastMousePos = event->pos();

	double sensitivity = 0.005;

	camera.orbitYaw(-delta.x() * sensitivity);
	camera.orbitPitch(-delta.y() * sensitivity);
	update();
}

void GraphWidget::wheelEvent(QWheelEvent *event)
{
	double numDegrees = event->angleDelta().y() / 120.0;

	if (2 <= dist*std::pow(0.9, numDegrees) && dist*std::pow(0.9, numDegrees) <= 50) {
		dist *= std::pow(0.9, numDegrees);
		camera.setDist(dist);
	}

	update();
}

GraphWidget::Polygon GraphWidget::pointToPixels(const R3Geometry::Point &p) const 
{	
	auto point = camera.projectPoint({p.X(), p.Y(), -p.Z()});
	GraphWidget::Polygon res; 
	res.vertex[0] = QPointF((0.5 + point.vertex[0].X())*w, (0.5 - point.vertex[0].Y())*h);
	res.n = 1;
	res.depth = point.depth;
	res.visible = point.visible;
	return res;
}

GraphWidget::Polygon GraphWidget::segToPixels(const R3Geometry::Point &p, const R3Geometry::Point &q) const 
{
	auto seg = camera.projectSegment({p.X(), p.Y(), -p.Z()}, {q.X(), q.Y(), -q.Z()});
	GraphWidget::Polygon res; 
	res.vertex[0] = QPointF((0.5 + seg.vertex[0].X())*w, (0.5 - seg.vertex[0].Y())*h);
	res.vertex[1] = QPointF((0.5 + seg.vertex[1].X())*w, (0.5 - seg.vertex[1].Y())*h);
	res.n = 2;
	res.depth = seg.depth;
	res.visible = seg.visible;
	return res;
}

GraphWidget::Polygon GraphWidget::triangleToPixels(const R3Geometry::Point &a, const R3Geometry::Point &b, const R3Geometry::Point &c) const 
{
	auto poly = camera.projectTriangle({a.X(), a.Y(), -a.Z()}, {b.X(), b.Y(), -b.Z()}, {c.X(), c.Y(), -c.Z()});
	GraphWidget::Polygon res; 
	for (int i = 0; i < poly.n; ++i)
		res.vertex[i] = QPointF((0.5 + poly.vertex[i].X())*w, (0.5 - poly.vertex[i].Y())*h);
	res.n = poly.n;
	res.depth = poly.depth;
	res.visible = poly.visible;
	res.color = converstion(light.compute(1./3*(a+b+c), (b-a).vect(c-a)));
	return res;
}

QColor GraphWidget::converstion(const PointLight::Color &color) const
{
	return QColor(std::round(color.X()), std::round(color.Y()), std::round(color.Z()));
}

QSize GraphWidget::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize GraphWidget::sizeHint () const
{
  return QSize (1000, 1000);
}
