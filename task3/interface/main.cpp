#include "utils/CommonDefs.h"

#include "interface/MainWindow.h"
#include <QApplication>
#include <sstream>
#include <iostream>

#define DEFAULT_A (-5.0)
#define DEFAULT_B (5.0)
#define DEFAULT_C (-5.0)
#define DEFAULT_D (5.0)
#define DEFAULT_NX (10)
#define DEFAULT_NY (10)
#define DEFAULT_MX (64)
#define DEFAULT_MY (64)
#define DEFAULT_K (0)

static
int parse_command_line (
		int argc, char **argv, 
		double &a, double &b, double &c, double &d, 
		int &nx, int &ny, int &mx, int &my, int &k)
{
    if (argc == 1)
        return 0;

    if (argc != 7)
    {
		std::cerr << "Wrong number of arguments!\n" << std::endl;
        return -1;
    }

	std::stringstream ss;
	for (int i = 1; i < argc; ++i)
		ss << argv[1] << " ";

    if (!(ss >> a >> b >> c >> d))
    {
        std::cerr << "Invalid value for bounds!\n" << std::endl;
        return -1;
	}

    if (!(ss >> nx >> ny))
    {
        std::cerr << "Invalid value for n!\n" << std::endl;
        return -1;
    }

	if (!(ss >> mx >> my))
	{
        std::cerr << "Invalid value for m!\n" << std::endl;
        return -1;
    }

 	if (!(ss >> k))
    {
        std::cerr << "Invalid value for k!\n" << std::endl;
        return -1;
    }


    double segment_length = std::min(b-a, d-c);
    if (segment_length < EPSILON_FOR_COMPARE)
    {
        std::cerr << "Too short segment: " << segment_length << std::endl;
        return -1;
    }

    if (nx <= 0 || ny <= 0 || mx <= 0 || my <= 0)
    {
		std::cerr << "Non-positive number of points!\n" << std::endl;
        return -1;
    }
    if (k <= 0)
    {
        std::cerr << "Non-positive number of function!\n" << std::endl;
        return -1;
    }

    return 0;
}

int main (int argc, char **argv) {
	double a = DEFAULT_A;
	double b = DEFAULT_B;
	double c = DEFAULT_C;
	double d = DEFAULT_D;

	int nx = DEFAULT_NX;
	int ny = DEFAULT_NY;
	int mx = DEFAULT_MX;
	int my = DEFAULT_MY;

	int k = DEFAULT_K;

	int err = parse_command_line (argc, argv, a, b, c, d, nx, ny, mx, my, k);
	if (err < 0)
		return -1;

	QApplication app (argc, argv);

	MainWindow w(a, b, c, d, nx, ny, mx, my, k);

	w.show();
	return app.exec();
}
