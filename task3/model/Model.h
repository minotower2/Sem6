#ifndef MODEL_H
#define MODEl_H

#include "gui/GraphWidget.h"
#include "gui/Data.h"
#include "utils/Geometry.h"
#include "resources/Storage.h"
#include "core/Approximation.h"

class Model : public QObject 
{
	Q_OBJECT
private:
	
	class Function {
		public:
		std::string name;
		double (*f)(double, double);
		Function(std::string, double (*)(double, double));
	};

	GraphWidget *graphWidget;
	
	int nx, ny;
	int mx, my;

	double ax, bx;
	double ay, by;

	std::vector<Function> functions;

	Storage *storage;
	Data *data;

	Approx::Approximator approximator;
public:
	Model(double, double, double, double, int, int, int, int, int, Storage *, Data *, QObject *);
	
	void errorUp(void);
	void errorDown(void);
	
	void changeFunc(void);
	void changeMode(void);
	void downTriangulation(void);
	void upTriangulation(void);

	void numberDown(void);
	void numberUp(void);
private:
	void updateBound(void);

	void remake(std::vector<R3Geometry::Point> &, double (Approx::Approximator::*) (double, double));
	
	void makeOrigin(std::vector<R3Geometry::Point> &);
	void makeApprox(std::vector<R3Geometry::Point> &);
	void makeResudial(std::vector<R3Geometry::Point> &);
	
	void update(void);
	void updateFull(void);
};

#endif // MODEl_H
