#ifndef GEOMETRY
#define GEOMETRY

#include <cmath>

#define EPS (1e-15)

namespace R2Geometry {

	class Point;
	class Vector;

	Vector operator* (const double c, const Vector &v);
	Point operator* (const double c, const Point &p);

	Vector RVector (const Point &p);
	Point RPoint (const Vector &v);
	
	class Point {
		double x;
		double y;

		public:

		Point (double xx = 0, double yy = 0) {
			if (fabs(xx) < EPS) x = 0;
			else x = xx;

			if (fabs(yy) < EPS) y = 0;
			else y = yy; 
		}
		
		Point (const Point &p) {
			x = p.x;
			y = p.y;
		}

		Point (Point &&p) {
			x = p.x; p.x = 0;
			y = p.y; p.y = 0;
		}

		Point & operator= (const Point &p) {
			x = p.x;
			y = p.y;
			return *this;
		}

		Point & operator= (Point &&p) {
			x = p.x; p.x = 0;
			y = p.y; p.y = 0;
			return *this;
		}
	
		double X (void) const {
			return x;
		}

		double Y (void) const {
			return y;
		}

		Vector operator- (const Point &) const;
		Point operator+ (const Vector &) const;
	};

	class Vector {
		double x;
		double y;

		public:

		Vector (double xx = 0, double yy = 0) {
			if (fabs(xx) < EPS) x = 0;
			else x = xx;

			if (fabs(yy) < EPS) y = 0;
			else y = yy;
		}

		Vector (const Vector &v) {
			x = v.x;
			y = v.y;
		}

		Vector (Vector &&v) {
			x = v.x; v.x = 0;
			y = v.y; v.y = 0;
		}

		Vector & operator= (const Vector &v) {
			x = v.x;
			y = v.y;
			return *this;
		}
	
		Vector & operator= (Vector &&v) {
			x = v.x; v.x = 0;
			y = v.y; v.y = 0;
			return *this;
		}

		double X (void) const{
			return x;
		}

		double Y (void) const{
			return y;
		}

		Vector operator+ (const Vector &v) const {
			return Vector(x+v.x, y+v.y);
		}

		Vector operator- (const Vector &v) const {
			return Vector(x-v.x, y-v.y);
		}

		Point operator+ (const Point &p) const {
			return Point(x+p.X(), y+p.Y());
		}

		Vector operator+ (void) const {
			return Vector(x, y);
		}

		Vector operator- (void) const {
			return Vector(-x, -y);
		}

		double operator* (const Vector v) const {
			return x*v.x+y*v.y;
		}

		Vector operator* (const double c) const {
			return Vector(c*x, c*y);
		}

		double mod (void) const{
			return sqrt(x*x);
		}

		Vector n (void) const {
			return Vector(-y, x);
		}

		Vector unit (void) const {
			return Vector(x/this->mod(), y/this->mod());
		}

		Vector & norm (void) {
			double len = this->mod();
			x /= len; y /= len;
			return *this;
		}
	};

}

namespace R3Geometry {

	class Point;
	class Vector;

	Vector operator* (const double c, const Vector &v);
	Point operator* (const double c, const Point &p);

	Vector RVector (const Point &p);
	
	Point RPoint (const Vector &v);
	
	class Point {
		double x;
		double y;
		double z;

		public:

		Point (double xx = 0, double yy = 0, double zz = 0) {
			if (fabs(xx) < EPS) x = 0;
			else x = xx;

			if (fabs(yy) < EPS) y = 0;
			else y = yy; 

			if (fabs(zz) < EPS) z = 0;
			else z = zz;
		}
		
		Point (const Point &p) {
			x = p.x;
			y = p.y;
			z = p.z;
		}

		Point (Point &&p) {
			x = p.x; p.x = 0;
			y = p.y; p.y = 0;
			z = p.z; p.z = 0;
		}

		Point & operator= (const Point &p) {
			x = p.x;
			y = p.y;
			z = p.z;
			return *this;
		}

		Point & operator= (Point &&p) {
			x = p.x; p.x = 0;
			y = p.y; p.y = 0;
			z = p.z; p.z = 0;
			return *this;
		}
	
		double X (void) const {
			return x;
		}

		double Y (void) const {
			return y;
		}

		double Z (void) const {
			return z;
		}

		Point operator+ (const Point &p) const
		{
			return Point(x+p.x, y+p.y, z+p.z);
		}

		Point operator+ (const Vector &) const;
		Vector operator- (const Point &) const;
		
		Point & operator+= (const Vector &v) {
			*this = *this + v;
			return *this;
		}
	
		double mod (void) const{
			return sqrt(x*x+y*y+z*z);
		}
	};

	class Vector {
		double x;
		double y;
		double z;

		public:

		Vector (double xx = 0, double yy = 0, double zz = 0) {
			if (fabs(xx) < EPS) x = 0;
			else x = xx;

			if (fabs(yy) < EPS) y = 0;
			else y = yy;

			if (fabs(zz) < EPS) z = 0;
			else z = zz;
		}

		Vector (const Vector &v) {
			x = v.x;
			y = v.y;
			z = v.z;
		}

		Vector (Vector &&v) {
			x = v.x; v.x = 0;
			y = v.y; v.y = 0;
			z = v.z; v.z = 0;
		}

		Vector & operator= (const Vector &v) {
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}
	
		Vector & operator= (Vector &&v) {
			x = v.x; v.x = 0;
			y = v.y; v.y = 0;
			z = v.z; v.z = 0;
			return *this;
		}

		double X (void) const{
			return x;
		}

		double Y (void) const{
			return y;
		}

		double Z (void) const{
			return z;
		}

		Vector operator+ (const Vector &v) const {
			return Vector(x+v.x, y+v.y, z+v.z);
		}

		Vector operator- (const Vector &v) const {
			return Vector(x-v.x, y-v.y, z-v.z);
		}

		Point operator+ (const Point &p) const {
			return Point(x+p.X(), y+p.Y(), z+p.Z());
		}

		Vector operator+ (void) const {
			return Vector(x, y, z);
		}

		Vector operator- (void) const {
			return Vector(-x, -y, -z);
		}

		double operator* (const Vector &v) const {
			return x*v.x+y*v.y+z*v.z;
		}

		Vector operator* (const double c) const {
			return Vector(c*x, c*y, c*z);
		}

		double mod (void) const{
			return sqrt(x*x+y*y+z*z);
		}

		Vector n (void) const {
			return Vector(-y, x, 0);
		}

		inline Vector unit (void) const {
			return Vector(x/this->mod(), y/this->mod(), z/this->mod());
		}

		inline Vector & norm (void) {
			double len = this->mod();
			x /= len; y /= len; z /= len;
			return *this;
		}

		inline Vector vect (const Vector &v) const {
			return Vector(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
		}
	};

}

class Quaternion {
	
	private:
	using Mat = struct Matrix3
	{
		double m[3][3];
		R3Geometry::Vector operator* (const R3Geometry::Vector &v) const noexcept
		{
			return {
				m[0][0]*v.X() + m[0][1]*v.Y() + m[0][2]*v.Z(),
				m[1][0]*v.X() + m[1][1]*v.Y() + m[1][2]*v.Z(),
				m[2][0]*v.X() + m[2][1]*v.Y() + m[2][2]*v.Z()
			};
		}
	};

	Mat m;
	public:
	double w, x, y, z;

	Quaternion (double w, double x, double y, double z);
	Quaternion (const R3Geometry::Vector &axis, double angle);

	void normalize ();

	Quaternion conjugate () const;
	Quaternion inverse () const;

	Quaternion operator* (const Quaternion &rhs) const;

	R3Geometry::Vector rotate (const R3Geometry::Vector &v) const;
};

#endif
