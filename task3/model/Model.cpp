#include "model/Model.h"

Model::Function::Function (std::string name, double(*f)(double, double)): name(name), f(f)
{}

Model::Model (double ax, double bx, double ay, double by, int nx, int ny, int mx, int my, int k, Storage *storage, Data *data, QObject *parent = nullptr) :
	QObject(parent),
	nx(nx), ny(ny),
	mx(mx), my(my),
	ax(ax), bx(bx),
	ay(ay), by(by),
	storage(storage), data(data),
	approximator(ax, bx, ay, by, nx, ny)
{	
	functions = std::vector<Function>
	{
		Function("f(x, y) = 1", [](double , double ) {return 1.0;}),
		Function("f(x, y) = x", [](double x, double ) {return x;}),
		Function("f(x, y) = y", [](double , double y) {return y;}),
		Function("f(x, y) = x+y", [](double x, double y) {return x+y;}),
		Function("f(x, y) = (x^2+y^2)^(1/2)", [](double x, double y) {return sqrt(x*x + y*y);}),
		Function("f(x, y) = x^2+y^2", [](double x, double y) {return x*x+ y*y;}),
		Function("f(x, y) = e^(x^2-y^2)", [](double x, double y) {return exp(x*x - y*y);}),
		Function("f(x, y) = 1/(25(x^2+y^2)+1)", [](double x, double y) {return 1/(25*(x*x+y*y)+1);})
	};
	
	data->k = k;
	data->name = functions[k].name;

	data->origin = true;
	data->piecePolynom = false;
	data->error = false;
	data->s = 0;
	data->nx = nx;
	data->ny = ny;
	data->mx = mx;
	data->my = my;
	data->p = 0;

	emit data->changed();

	approximator.setFunction(functions[data->k].f);
	updateBound();

	updateFull();
}

void Model::changeFunc(void)
{
	data->k = (data->k+1)%functions.size();
	data->name = functions[data->k].name;

	data->p = 0;
	
	emit data->changed();

	approximator.setFunction(functions[data->k].f);

	updateBound();

	update();
}

void Model::updateBound(void)
{
	std::pair<double, double> buff;

	buff = approximator.minMaxChangeable(&Approx::Approximator::origin);
	data->originAbsMax = std::max(fabs(buff.first), fabs(buff.second));
	
	if (data->piecePolynom) {
		buff = approximator.minMaxChangeable(&Approx::Approximator::approxPiecePolynom);
		data->piecePolynomAbsMax = std::max(fabs(buff.first), fabs(buff.second));
	}	
	if (data->error) {
		buff = approximator.minMaxChangeable(&Approx::Approximator::errorPiecePolynom);
		data->errorPiecePolynomAbsMax = std::max(fabs(buff.first), fabs(buff.second));
	}

	emit data->changed();
}

void Model::changeMode(void)
{
	if(data->origin) {
		data->origin = !(data->piecePolynom = data->origin);
	}	
	else if(data->piecePolynom) {
		data->piecePolynom = !(data->error = data->piecePolynom);
	}
	else if(data->error) {
		data->error = !(data->origin = data->error);
	}		

	emit data->changed();

	updateBound();

	update();
}

void Model::downTriangulation(void)
{
	if (mx*my >= 16) {
		mx /= 2;
		my /= 2;

		data->mx = mx;
		data->my = my;

		emit data->changed();
	
		updateFull();
	}
}

void Model::upTriangulation(void)
{
	if (mx*my <= 128*128) {
		mx *= 2;
		my *= 2;
	
		data->mx = mx;
		data->my = my;

		emit data->changed();

		updateFull();
	}
}

void Model::errorUp()
{
	if (data->p < 10) {
		++(data->p);
		approximator.setError(data->p*0.1*data->originAbsMax);
		
		emit data->changed();

		updateBound();
		update();
	}

}

void Model::errorDown()
{
	if (data->p > -10) {
		--(data->p);
		approximator.setError(data->p*0.1*data->originAbsMax);

		emit data->changed();

		updateBound();
		update();
	}
}

void Model::numberDown(void)
{
	if (nx*ny >= 100) {
		nx /= 2;
		ny /= 2;

		data->nx = nx;
		data->ny = ny;

		emit data->changed();
	
		approximator.setN(nx, ny);

		updateBound();
		update();
	}
}

void Model::numberUp(void)
{
	if (nx*ny <= 25*100000 && nx*ny <= 25*25) {
		nx += 1;
		ny += 1;

		data->nx = nx;
		data->ny = ny;

		emit data->changed();
	
		approximator.setN(nx, ny);

		updateBound();
		update();
	}

}

void Model::remake(std::vector<R3Geometry::Point> &points, double (Approx::Approximator::*f) (double, double))
{
	double dx = (bx-ax)/(mx-1), dy = (by-ay)/(my-1);
	for (int i = 0; i < mx; ++i) {
		for (int j = 0; j < my; ++j) {
			points[i*my + j] = R3Geometry::Point(ax+i*dx, ay+j*dy, (approximator.*f)(ax+i*dx, ay+j*dy) );
		}
	}
}

void Model::makeOrigin(std::vector<R3Geometry::Point> &points)
{
	remake(points, &Approx::Approximator::origin);
}

void Model::makeApprox(std::vector<R3Geometry::Point> &points)
{
	remake(points, &Approx::Approximator::approxPiecePolynom);
}

void Model::makeResudial(std::vector<R3Geometry::Point> &points)
{
	remake(points, &Approx::Approximator::errorPiecePolynom);
}

void Model::update(void)
{
	auto mesh = storage->getAccess();

	if (data->origin)
		makeOrigin(mesh->points);
	if (data->piecePolynom)
		makeApprox(mesh->points);
	if (data->error)
		makeResudial(mesh->points);
	storage->update();
}

void Model::updateFull(void)
{
	auto mesh = storage->getAccess();

	mesh->points.resize(mx*my);

	if (data->origin)
		makeOrigin(mesh->points);
	if (data->piecePolynom)
		makeApprox(mesh->points);
	if (data->error)
		makeResudial(mesh->points);

	mesh->triangles.resize((mx-1)*(my-1)*2);
	for (int i = 0; i < mx-1; ++i) {
		for (int j = 0; j < my-1; ++j) {
			mesh->triangles[2*(i*(my-1) + j)] = Storage::Triangle(i*my+j, i*my+j+1, (i+1)*my+j, 0);
			mesh->triangles[2*(i*(my-1) + j)+1] = Storage::Triangle(i*my+j+1, (i+1)*my+j+1, (i+1)*my+j, 0);
		}
	}

	storage->update();
}
