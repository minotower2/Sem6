#include "scene3D.h"
#include <QApplication>
// #include <QMainWindow> // Не используется напрямую
// #include <QMenuBar>    // Не используется
// #include <QMessageBox> // Не используется
// #include <QPainter>    // Не используется
// #include <QVBoxLayout> // Не используется

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Scene3D scene1(argv);
    scene1.setWindowTitle("Аппроксимация сплайнами"); // Изменен заголовок
    if (scene1.parse_command_line(argc, argv)) { // Этот вызов все равно будет избыточным, если конструктор его уже сделал успешно
        qWarning("Wrong input arguments or error during re-parsing!"); // Сообщение изменено для ясности
        return -1;
    }

    scene1.resize(500, 500); // Уменьшим для удобства, было 5000x5000
    scene1.move(0,0);
    scene1.show();
    return app.exec();
}
