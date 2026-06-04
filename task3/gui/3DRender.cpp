#include "3DRender.h"

Camera::Camera(R3Geometry::Point target, Quaternion orien, double f, double d, double asp): 
	target(target), orientation(orien), focus(f), dist(d),
	screenHalfHeight(1/sqrt(3)), aspect(asp),
	near(f), far(100)
{
	orientation.normalize();
	W = near/(screenHalfHeight*aspect);
	H = near/(screenHalfHeight);
	max = screenHalfHeight/(2*near);
}

bool Camera::check(const R3Geometry::Vector &v) const
{
	return v.Z() <= -near && v.Z() >= - far && fabs(v.X()) <= -v.Z()*aspect*max && fabs(v.Y()) <= -v.Z()*max;
}

bool Camera::intersect(R3Geometry::Vector &va, R3Geometry::Vector &vb) const
{
	double t0 = 0, t1 = 1, t;
	R3Geometry::Vector dir = vb - va;

	double fa = 0, fb = 0;

	for (int i = 0; i < 6; ++i) {
		switch(i) {	
			case 0:
			fa = va.Z() + near;
		   	fb = vb.Z() + near;
			break;
			
			case 1:
			fa = -va.Z() - far;
		   	fb = -vb.Z() - far;
			break;

			case 2:
			fa = va.X() + va.Z()*aspect*max;
		   	fb = vb.X() + vb.Z()*aspect*max;
			break;

			case 3:
			fa = -va.X() + va.Z()*aspect*max;
		   	fb = -vb.X() + vb.Z()*aspect*max;
			break;

			case 4:
			fa = va.Y() + va.Z()*max;
		   	fb = vb.Y() + vb.Z()*max;
			break;

			case 5:
			fa = -va.Y() + va.Z()*max;
		   	fb = -vb.Y() + vb.Z()*max;
			break;
		}

		if (fa > 0 && fb > 0)
			return false;
				
		if (fa <= 0 && fb <= 0)
			continue;

		t = -fa / (fb - fa);

		if (fb-fa < 0.0)
			t0 = std::max(t0, t);
		else
			t1 = std::min(t1, t);

		if (t0 > t1)
			return false;
	}
	
	auto v = va + t0*dir;
	vb = va + t1*dir;
	va = v;

	return true;
}



R3Geometry::Vector Camera::toView(const R3Geometry::Point &x) const 
{
	return orientation.conjugate().rotate(x - position());
}

R2Geometry::Point Camera::toProj(const R3Geometry::Vector &v) const 
{
	double div = 1/v.Z(); 
	double x = (-W) * v.X()*div;
	double y = (-H) * v.Y()*div;
	return R2Geometry::Point(x, y);
}

Camera::Polygon Camera::projectPoint(const R3Geometry::Point &p) const
{
	auto vp = toView(p);
	Camera::Polygon res;
	res.vertex[0] = {0, 0};
	res.n = 1;
	res.depth = p.mod();
	res.visible = check(vp);

	if (res.visible)
		res.vertex[0] = toProj(vp);

	return res;
}

Camera::Polygon Camera::projectSegment(const R3Geometry::Point &a, const R3Geometry::Point &b) const
{
	auto va = toView(a);
	auto vb = toView(b);

	Camera::Polygon res;
	res.vertex[0] = {0, 0};
	res.vertex[1] = {0, 0};
	res.n = 2;

	res.depth = (1./2 *(a+b)).mod();

	res.visible = false;

	if (!intersect(va, vb))
		return res;

	res.vertex[0] = toProj(va);
	res.vertex[1] = toProj(vb);
	res.visible = true;

//	std::cerr << "a: " << res.a.X() << " " << res.a.Y() << "\nb: " << res.b.X() << " " << res.b.Y() << std::endl; 

	return res;
}

Camera::Polygon Camera::projectTriangle(const R3Geometry::Point &a, const R3Geometry::Point &b,  const R3Geometry::Point &c) const
{
	auto va = toView(a), vb = toView(b), vc = toView(c);

	Camera::Polygon res;

	res.vertex[0] = {0, 0};
	res.vertex[1] = {0, 0};
	res.vertex[2] = {0, 0};

	res.n = 3;
	res.visible = false;

	res.depth = (1./3 *(a+b+c)).mod();

	bool aIn = check(va), bIn = check(vb), cIn = check(vc);
	if (!aIn && !bIn && !cIn)
		return res;

	res.visible = true;
	if (aIn && bIn && cIn) {	
		res.vertex[0] = toProj(va);
		res.vertex[1] = toProj(vb);
		res.vertex[2] = toProj(vc);
	}
	else if (!aIn && bIn && cIn) {
		res.n = 4;
		auto tmp = va;
		intersect(tmp, vb);
		intersect(va, vc);
		res.vertex[0] = toProj(tmp);
		res.vertex[1] = toProj(vb);
		res.vertex[2] = toProj(vc);
		res.vertex[3] = toProj(va);
	}
	else if (aIn && !bIn && cIn) {
		res.n = 4;
		auto tmp = vb;
		intersect(tmp, va);
		intersect(vb, vc);
		res.vertex[0] = toProj(tmp);
		res.vertex[1] = toProj(va);
		res.vertex[2] = toProj(vc);
		res.vertex[3] = toProj(vb);

	}
	else if (aIn && bIn && !cIn) {
		res.n = 4;
		auto tmp = vc;
		intersect(tmp, va);
		intersect(vc, vb);
		res.vertex[0] = toProj(tmp);
		res.vertex[1] = toProj(va);
		res.vertex[2] = toProj(vb);
		res.vertex[3] = toProj(vc);
	}
	
	else if (!aIn && !bIn && cIn) {
		auto tmp = vc;
		intersect(va, tmp);
		tmp = vc;
		intersect(vb, tmp);

		res.vertex[0] = toProj(va);
		res.vertex[1] = toProj(vb);
		res.vertex[2] = toProj(vc);
	}
	
	else if (!aIn && bIn && !cIn) {
		auto tmp = vb;
		intersect(va, tmp);
		tmp = vb;
		intersect(vc, tmp);

		res.vertex[0] = toProj(va);
		res.vertex[1] = toProj(vb);
		res.vertex[2] = toProj(vc);

	}	
	else if (aIn && !bIn && !cIn) {
		auto tmp = va;
		intersect(vc, tmp);
		tmp = va;
		intersect(vb, tmp);

		res.vertex[0] = toProj(va);
		res.vertex[1] = toProj(vb);
		res.vertex[2] = toProj(vc);
	}
	
	return res;
}

void Camera::orbitYaw(double angle)
{
	orbitGlobal(R3Geometry::Vector(0,0,-1), angle);
}

void Camera::orbitPitch(double angle)
{
	orbitLocal(R3Geometry::Vector(1,0,0), angle);
}

void Camera::moveForward(double dist)
{
	move({0, 0, -1}, dist);
}

void Camera::moveRight(double dist)
{
	move({1, 0, 0}, dist);
}

void Camera::moveUp(double dist)
{
	move({0, 1, 0}, dist);
}



void Camera::move(const R3Geometry::Vector &localDir, double dist)
{
	target += dist * orientation.rotate(localDir.unit());
}

void Camera::orbitGlobal(const R3Geometry::Vector &axis, double phi)
{
	Quaternion q(axis.unit(), phi);
	
	orientation = q * orientation;
	orientation.normalize();
}

void Camera::orbitLocal(const R3Geometry::Vector &axis, double phi)
{
	Quaternion q(axis.unit(), phi);
	
	orientation = orientation * q;
	orientation.normalize();
}


void Camera::setDist(double d) 
{
	dist = d;
}

void Camera::setAspect(double asp)
{
	aspect = asp;
	W = near/(screenHalfHeight*aspect);
}

R3Geometry::Point Camera::position(void) const
{
	return target +  dist*orientation.rotate({0,0,1});
}

std::string Camera::orien (void) const 
{
	return 	"(" + std::to_string(orientation.w) + ", "+ std::to_string(orientation.x) + ", "+ std::to_string(orientation.y) + ", "+ std::to_string(orientation.z) + ")\n";
}


PointLight::Color PointLight::compute(const R3Geometry::Point &point, const R3Geometry::Vector &normal) const
{
	R3Geometry::Vector L(point-position);
	double dist = L.mod();
    double attenuation = 1.0 / (1.0 + 0.09 * dist + 0.032 * dist * dist);
	return (ambient + intensity * std::max(0.0, (normal.unit() * (1/dist * L))) * attenuation)*color;
}
