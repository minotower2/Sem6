#ifndef STORAGE_H
#define STORAGE_H

#include <vector>
#include <QObject> 
#include <QPointF>
#include <QColor>

#include "utils/Geometry.h"

class Storage : public QObject
{
	Q_OBJECT
	public:

	using Polygon = struct 
	{
	    QPointF vertex[4];
		int n;
		double depth;
		QColor color;
	    bool visible;
	};

	class Triangle
	{
		public:
		int a, b, c;
		double depth;
		
		Triangle (int a = 0, int b = 0, int c = 0, double depth = 0);
		Triangle (const Storage::Triangle &tr);
		Triangle (Storage::Triangle && tr);
		Triangle & operator= (const Triangle &tr);
		Triangle & operator= (Triangle && tr);

	};

	private:
	class Mesh
	{
		public:
		
		std::vector<R3Geometry::Point> points;
		std::vector<Triangle> triangles; 
		Mesh(void)
		{
		}

		Mesh& operator= (const Mesh& mesh)
		{
			points = mesh.points;
			triangles = mesh.triangles;
			return *this;
		}
	};

	bool active;
	Mesh currentFrame;
	Mesh readyFrame;
	Mesh newFrame;

	std::vector<Polygon> polygons;
	public:
	Storage(QObject *object = nullptr);

	Mesh * renderrAccess(void);
	Mesh * getAccess(void);
	std::vector<Polygon> & poly(void);
	void update(void);

	signals:
	void updated();
};

#endif //STORAGE_H
