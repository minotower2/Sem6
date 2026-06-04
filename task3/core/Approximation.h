#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include "utils/Geometry.h"
#include "utils/CommonDefs.h"
#include <algorithm>
#include <cassert>
#include <vector>

#define TIMER (0)

namespace Approx {

	class Approximator {
		double ax, bx;
		double ay, by;
		int nx, ny;
		double p;

		double (*f)(double, double);

		double pixel;
		R2Geometry::Point mid;

		std::vector<double> points_x;
		std::vector<double> points_y;

		std::vector<double> values;
		std::vector<double> coeff;

		void makeApproxPiecePolynom(void);

		public:
		Approximator(double, double, double, double, int, int, double (*)(double, double) = [](double, double) {return 0.0;});
	
		double origin(double x, double y);
		double approxPiecePolynom(double x, double y);
		double errorPiecePolynom(double x, double y);

		std::pair<double, double> minMaxChangeable(double (Approximator::*)(double, double));

		std::vector<std::pair<double, double>> minMax(void);

		void setFunction(double (*)(double, double));
		void setError(double);
		void setN(int, int);
		void setPixel(double);
	};
}

#endif //APPROXIMATION_H
