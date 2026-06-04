#include "resources/Storage.h"

Storage::Triangle::Triangle (int a, int b, int c, double depth): a(a), b(b), c(c), depth(depth)
{}

Storage::Triangle::Triangle (const Storage::Triangle &tr)
{
	a = tr.a; b = tr.b; c = tr.c; depth = tr.depth;
}
Storage::Triangle::Triangle (Storage::Triangle && tr)
{
	a = tr.a; b = tr.b; c = tr.c; depth = tr.depth;
	tr.a = 0; tr.b = 0; tr.c = 0; tr.depth = 0;
}

Storage::Triangle & Storage::Triangle::operator= (const Storage::Triangle &tr)
{
	a = tr.a; b = tr.b; c = tr.c; depth = tr.depth;
	return *this;
}

Storage::Triangle & Storage::Triangle::operator= (Storage::Triangle && tr)
{
	a = tr.a; b = tr.b; c = tr.c; depth = tr.depth;
	tr.a = 0; tr.b = 0; tr.c = 0; tr.depth = 0;
	return *this;
}

Storage::Storage(QObject *object):
	QObject(object)
{
}

Storage::Mesh * Storage::renderrAccess(void)
{
	if (active)
		currentFrame = readyFrame;

	return &currentFrame;
}

Storage::Mesh * Storage::getAccess(void)
{
	active = false;
	return &newFrame;
}

void Storage::update(void) 
{
	readyFrame = newFrame;
	active = true;
	emit updated();
}

std::vector<Storage::Polygon> & Storage::poly(void)
{
	return polygons;
}
