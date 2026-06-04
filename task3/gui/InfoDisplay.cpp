#include "gui/InfoDisplay.h"

InfoDisplay::InfoDisplay (QWidget *parent, Data *data):
	QWidget(parent),
	data(data),
	bgColor(Qt::white)
{
	connect(data, &Data::changed, this, &InfoDisplay::dataChanged);

	table = new QVBoxLayout(this);

	label = new QLabel(this);
	label->setAlignment(Qt::AlignTop | Qt::AlignLeft);
//	label->setTextInteractionFlags(Qt::TextSelectableByMouse);

	table->addWidget(label);
	table->addStretch();
}


void InfoDisplay::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);

	int w = width();
	int h = height();
	if (w == 0 || h == 0)
		return;

	painter.setBrush(QBrush(bgColor));
	painter.drawRect(0, 0, w, h);
	painter.setBrush(Qt::NoBrush);
	
	painter.setRenderHint(QPainter::Antialiasing, true);
	drawInformation();
	drawConsole();
	event->accept();
}

void InfoDisplay::drawConsole(void)
{
	static bool flag = false;
	static int lines = 0;
	if (flag) {

		std::cout << "\033[" << lines << "A";
		lines = 0;
		std::cout << "\033[2K\r" << "k = " << data->k << " " + data->name + "\n"; 
		++lines;
		if (data->origin) {
			std::cout << "\033[2K\r" << "max |f_origin(x)| = " << data->originAbsMax << "\n";
			++lines;
		}
		if (data->piecePolynom) {
			std::cout << "\033[2K\r" << "max |f_piece_polynom(x)| = " << data->piecePolynomAbsMax << "\n";
			++lines;
		}
		if (data->error) {
			std::cout << "\033[2K\r" << "max |f_piece_polynom_error(x)| = " << data->errorPiecePolynomAbsMax << "\n";
			++lines;
		}

		std::cout.flush();
	} else {
		std::cout << "k = " << data->k << " " + data->name + "\n"; 
		++lines;
		if (data->origin) {
			std::cout << "max |f_origin(x)| = " << data->originAbsMax << "\n";
			++lines;
		}
		if (data->piecePolynom) {
			std::cout << "max |f_piece_polynom(x)| = " << data->piecePolynomAbsMax << "\n";
			++lines;
		}
		if (data->error) {
			std::cout << "max |f_piece_polynom_error(x)| = " << data->errorPiecePolynomAbsMax << "\n";
			++lines;
		}

		flag = true;
	}

}
void InfoDisplay::drawInformation(void) 
{

	QString text;

	text += QString("k = %1 %2\n")
		.arg(data->k)
		.arg(QString::fromStdString(data->name));

	if (data->origin)
		text += QString("max |f_origin(x)|:\n %1\n").arg(data->originAbsMax);

	if (data->piecePolynom)
		text += QString("max |f_piece_polynom(x)|:\n %1\n").arg(data->piecePolynomAbsMax);

	if (data->error) {
		text += QString("max |f_piece_polynom_error(x)|:\n %1\n").arg(data->errorPiecePolynomAbsMax);
	}
	
	text += QString("\ns = %1\nn_x = %2    n_y = %3\nm_x = %4   m_y = %5\np = %6")
		.arg(data->s)
		.arg(data->nx)
		.arg(data->ny)
		.arg(data->mx)
		.arg(data->my)
		.arg(data->p);

	label->setText(text);
}
