#include "Geometry.h"

namespace R2Geometry {
	Vector operator*(const double c, const Vector& v) 
	{
		return Vector(c*v.X(), c*v.Y());
	}

	Point operator* (const double c, const Point& p) 
	{
		return Point(c*p.X(), c*p.Y());
	}

	Vector RVector (const Point& p) 
	{
		return Vector(p.X(), p.Y());
	}

	Point RPoint (const Vector& v) 
	{
		return Point(v.X(), v.Y());
	}

	Vector Point::operator- (const Point &p) const
	{
		return Vector(x-p.x, y-p.y);
	}

	Point Point::operator+ (const Vector &v) const
	{
		return Point(x+v.X(), y+v.Y());
	}
	
}

namespace R3Geometry {

	Vector operator*(const double c, const Vector& v) {
		return Vector(c*v.X(), c*v.Y(), c*v.Z());
	}

	Point operator* (const double c, const Point& p) {
		return Point(c*p.X(), c*p.Y(), c*p.Z());
	}

	Vector RVector (const Point& p) {
		return Vector(p.X(), p.Y(), p.Z());
	}

	Point RPoint (const Vector& v) {
		return Point(v.X(), v.Y(), v.Z());
	}
	
	Point Point::operator+ (const Vector &v) const
	{
		return Point(x+v.X(), y+v.Y(), z+v.Z());
	}

	Vector Point::operator- (const Point &p) const
	{
		return Vector(x-p.x, y-p.y, z-p.z);
	}
	
}

Quaternion::Quaternion(double www, double xxx, double yyy, double zzz):
	w(www), x(xxx), y(yyy), z(zzz)
{
	double mod = sqrt(w*w+x*x+y*y+z*z);
	w /= mod;
	x /= mod;
	y /= mod;
	z /= mod;

	const double xx = x*x;
	const double yy = y*y;
	const double zz = z*z;

	const double xy = x*y;
	const double xz = x*z;
	const double yz = y*z;

	const double wx = w*x;
	const double wy = w*y;
	const double wz = w*z;
		
	m.m[0][0] = 1.0 - 2.0*(yy + zz);
	m.m[0][1] = 2.0*(xy - wz);
	m.m[0][2] = 2.0*(xz + wy);

	m.m[1][0] = 2.0*(xy + wz);
	m.m[1][1] = 1.0 - 2.0*(xx + zz);
	m.m[1][2] = 2.0*(yz - wx);

	m.m[2][0] = 2.0*(xz - wy);
	m.m[2][1] = 2.0*(yz + wx);
	m.m[2][2] = 1.0 - 2.0*(xx + yy);

	/*
	if (!(mod > 1e-15)) {
		w /= mod;
		x /= mod;
		y /= mod;
		z /= mod;

		const double xx = x*x;
		const double yy = y*y;
		const double zz = z*z;

		const double xy = x*y;
		const double xz = x*z;
		const double yz = y*z;

		const double wx = w*x;
		const double wy = w*y;
		const double wz = w*z;
		
		m.m[0][0] = 1.0 - 2.0*(yy + zz);
		m.m[0][1] = 2.0*(xy - wz);
		m.m[0][2] = 2.0*(xz + wy);

		m.m[1][0] = 2.0*(xy + wz);
		m.m[1][1] = 1.0 - 2.0*(xx + zz);
		m.m[1][2] = 2.0*(yz - wx);

		m.m[2][0] = 2.0*(xz - wy);
		m.m[2][1] = 2.0*(yz + wx);
		m.m[2][2] = 1.0 - 2.0*(xx + yy);
	}	
	*/
}

Quaternion::Quaternion(const R3Geometry::Vector& axis, double angle)
{
	R3Geometry::Vector n = axis.unit();

	double half = angle * 0.5;
	double s = sin(half);

	w = cos(half);
	x = n.X() * s;
	y = n.Y() * s;
	z = n.Z() * s;

	const double xx = x*x;
	const double yy = y*y;
	const double zz = z*z;

	const double xy = x*y;
	const double xz = x*z;
	const double yz = y*z;

	const double wx = w*x;
	const double wy = w*y;
	const double wz = w*z;
		
	m.m[0][0] = 1.0 - 2.0*(yy + zz);
	m.m[0][1] = 2.0*(xy - wz);
	m.m[0][2] = 2.0*(xz + wy);

	m.m[1][0] = 2.0*(xy + wz);
	m.m[1][1] = 1.0 - 2.0*(xx + zz);
	m.m[1][2] = 2.0*(yz - wx);

	m.m[2][0] = 2.0*(xz - wy);
	m.m[2][1] = 2.0*(yz + wx);
	m.m[2][2] = 1.0 - 2.0*(xx + yy);
}

void Quaternion::normalize()
{
	double mod = sqrt(w*w+x*x+y*y+z*z);
	if (mod < 1e-15)
		return;
	w /= mod;
	x /= mod;
   	y /= mod;
   	z /= mod;
}

Quaternion Quaternion::conjugate() const
{
	return Quaternion(w, -x, -y, -z);
}

Quaternion Quaternion::inverse() const
{
	double mod2 = w*w+x*x+y*y+z*z;
	return Quaternion(w/mod2, -x/mod2, -y/mod2, -z/mod2);
}

Quaternion Quaternion::operator*(const Quaternion& rhs) const
{
	return Quaternion(
			w*rhs.w - (x*rhs.x+y*rhs.y+z*rhs.z), 
			w*rhs.x +rhs.w*x + (y*rhs.z - z*rhs.y),
			w*rhs.y +rhs.w*y + (z*rhs.x - x*rhs.z),
			w*rhs.z +rhs.w*z + (x*rhs.y - y*rhs.x)
			);
}

R3Geometry::Vector Quaternion::rotate(const R3Geometry::Vector& v) const
{
	return m*v;
	/*
	R3Geometry::Vector qv(x, y, z);
	
	R3Geometry::Vector t = 2.0 * qv.vect(v);
	return v + w*t + qv.vect(t);
	*/
}

