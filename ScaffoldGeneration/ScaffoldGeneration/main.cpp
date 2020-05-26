#include "scaffoldgeneration.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ScaffoldGeneration w;
	w.setGeometry(100, 100, 640, 480);
	w.show();
	return a.exec();
}
