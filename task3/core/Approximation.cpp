#include "Approximation.h"
#include <sys/time.h>
#include <iostream>

namespace Approx {

	Approximator::Approximator(double ax, double bx, double ay, double by, int nx, int ny, double (*f)(double, double) ): 
		ax(ax), bx(bx), 
		ay(ay), by(by),
		nx(nx), ny(ny), 
		p(0), f(f), 
		pixel(0.01), mid(R2Geometry::Point(ax+nx*(bx-ax)/(2*(nx-1)), ay+ny*(by-ay)/(2*(ny-1)))),
		points_x(std::vector<double>(nx)), points_y(std::vector<double>(ny)),
		values(std::vector<double>(nx*ny)), coeff(std::vector<double>(nx*ny))
	{
		for (int i = 0; i < nx; ++i)
			points_x[i] = (ax+bx)/2 + ((bx-ax)/2) * cos(M_PI*(i+0.5)/nx);
		for (int j = 0; j < ny; ++j)
			points_y[j] = (ay+by)/2 + ((by-ay)/2) * cos(M_PI*(j+0.5)/ny);

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				values[i*ny+j] = f(points_x[i], points_y[j]) + (i == nx/2 && j == ny/2? p : 0);
			}
		}

		mid = R2Geometry::Point(points_x[nx/2], points_y[ny/2]);

		makeApproxPiecePolynom();
	}

	void Approximator::makeApproxPiecePolynom(void)
	{
		

		coeff = values;

		for (int j = 0; j < ny; j++) {
			for (int i = 1; i < nx; ++i) {
				for (int k = nx-1; k >= i; --k) {
					coeff[k*ny+j] = (coeff[k*ny+j] - coeff[(k-1)*ny+j]) / (points_x[k] - points_x[k - i]);
				}
			}
		}

		for (int i = 0; i < nx; ++i) {
			for (int j = 1; j < ny; ++j) {
				for (int k = ny-1; k >= j; --k) {
					coeff[i*ny+k] = (coeff[i*ny+k] - coeff[i*ny+(k-1)]) / (points_y[k] - points_y[k - j]);
				}
			}
		}


		return;
	}	
	
	double Approximator::origin(double x, double y) 
	{
		return f(x, y) + (fabs(x-mid.X()) < pixel && fabs(y-mid.Y()) < pixel ? p : 0);
	}

	double Approximator::approxPiecePolynom(double x, double y)
	{
		double res = 0, buff = 0;

		buff = coeff[(nx-1)*ny + ny-1];
		for (int j = ny-2; j >= 0; --j) {
			buff = buff*(y - points_y[j])+ coeff[(nx-1)*ny + j];
		}
		res = buff;

	 	for (int i = nx-2; i >= 0; --i) {
			buff = coeff[i*ny + ny-1];
			for (int j = ny-2; j >= 0; --j) {
				buff = buff*(y - points_y[j])+ coeff[i*ny + j];
			}

			res = res*(x-points_x[i]) + buff;
		}

		return res;
	}

	double Approximator::errorPiecePolynom(double x, double y)
	{
		return fabs(f(x, y)-approxPiecePolynom(x, y));
	}

	void Approximator::setFunction(double (*func)(double, double))
	{
		f = func;
		p = 0;
		
		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				values[i*ny+j] = f(points_x[i], points_y[j]) + (i == nx/2 && j == ny/2? p : 0);
			}
		}

		makeApproxPiecePolynom();
	}

	void Approximator::setError(double pp)
	{
		p = pp;
		values[(nx/2) * ny + (ny/2)] = f(points_x[nx/2], points_y[ny/2]) + p;

		makeApproxPiecePolynom();
	}

	
	void Approximator::setN(int nxx, int nyy)
	{
		nx = nxx;
		ny = nyy;
	
		points_x.resize(nx);
		points_y.resize(ny);

		values.resize(nx*ny);
		coeff.resize(nx*ny);

		for (int i = 0; i < nx; ++i)
			points_x[i] = (ax+bx)/2 + ((bx-ax)/2) * cos(M_PI*(i+0.5)/nx);
		for (int j = 0; j < ny; ++j)
			points_y[j] = (ay+by)/2 + ((by-ay)/2) * cos(M_PI*(j+0.5)/ny);

		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				values[i*ny+j] = f(points_x[i], points_y[j]) + (i == nx/2 && j == ny/2? p : 0);
			}
		}

		mid = R2Geometry::Point(points_x[nx/2], points_y[ny/2]);

		makeApproxPiecePolynom();
	}

	void Approximator::setPixel(double Pixel) 
	{
		pixel = Pixel;
	}

	std::pair<double, double> Approximator::minMaxChangeable(double (Approximator::*ff)(double, double))
	{	
		std::pair<double, double> minMax;
		minMax.first = (this->*ff)(ax, ay);
		minMax.second = (this->*ff)(ax, ay);
		double xDelta = 0.01;
		double yDelta = 0.01;
		for (double x = ax; x - bx < EPSILON_FOR_COMPARE; x += xDelta)
		{
			for (double y = ay; y - by < EPSILON_FOR_COMPARE; y += yDelta) {
				double z = (this->*ff)(x, y);
				if (z < minMax.first)
					minMax.first = z;
				if (z > minMax.second)
					minMax.second = z;
			}
		}
		return minMax;
	}

	std::vector<std::pair<double, double>> Approximator::minMax(void)
	{	
		std::vector<std::pair<double, double>> minMax{};
		minMax.push_back(minMaxChangeable(&Approximator::origin));
		minMax.push_back(minMaxChangeable(&Approximator::approxPiecePolynom));
		minMax.push_back(minMaxChangeable(&Approximator::errorPiecePolynom));
		return minMax;
	}
}
