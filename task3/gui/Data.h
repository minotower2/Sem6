#ifndef DATA_H
#define DATA_H

#include <string>

#include <QObject> 

class Data : public QObject
{
	Q_OBJECT
public:
	int k;
	std::string name;
	double originAbsMax;
	double piecePolynomAbsMax;
	double errorPiecePolynomAbsMax;

	bool origin, piecePolynom, error;
	int s, nx, ny, mx, my, p;

	Data(QObject *parent = nullptr);
signals:
	void changed();
};

#endif // DATA_H
