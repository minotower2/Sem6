#ifndef RENDERER_H
#define RENDERER_H

#include "utils/Geometry.h"
#include <iostream>
#include <cmath>
#include <utility>

class Camera 
{
	using Polygon = struct 
	{
		R2Geometry::Point vertex[4];
		int n;
		double depth;
	    bool visible;
	};

	R3Geometry::Point target;
	
	Quaternion orientation;
	double focus;
	double dist;

	double screenHalfHeight;
	double aspect;

	double W, H;

	double near, far;
	double max;
	public:

	Camera(R3Geometry::Point, Quaternion, double, double, double);

	Polygon projectPoint(const R3Geometry::Point &) const;
	Polygon projectSegment(const R3Geometry::Point &, const R3Geometry::Point &) const;
	Polygon projectTriangle(const R3Geometry::Point &, const R3Geometry::Point &,  const R3Geometry::Point &) const;

	void orbitYaw(double);
	void orbitPitch(double);

	void moveForward(double);
	void moveRight(double);
	void moveUp(double);
	
	void setDist(double);
	void setAspect(double);

	std::string orien (void) const;

	private:

	bool check(const R3Geometry::Vector &) const;
	bool intersect(R3Geometry::Vector &, R3Geometry::Vector &) const;

	R3Geometry::Point position() const;

	R3Geometry::Vector toView(const R3Geometry::Point &) const;
	R2Geometry::Point toProj(const R3Geometry::Vector &) const;

	void orbitGlobal(const R3Geometry::Vector &, double);
	void orbitLocal(const R3Geometry::Vector &, double);

	void move(const R3Geometry::Vector &, double);
};

class Light {
	public:

	using Color = R3Geometry::Point;

	Color color = {255,0,0};
	double intensity = 1.0;
	
	virtual Color compute(const R3Geometry::Point &, const R3Geometry::Vector &) const = 0;
	virtual ~Light() = default;
};

class PointLight : public Light 
{
	public:
    double ambient = 0.1;
	R3Geometry::Point position;

	PointLight (const R3Geometry::Point &position): position(position) {}
	Color compute(const R3Geometry::Point &, const R3Geometry::Vector &) const override;
};


#endif //3DRENDER_H
