#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGroupBox>
#include "sundial.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    QVector<double> dt;
    QVector<double> day;
    QVector<double> dt_between;
    QVector<double> dt_summer;
    QVector<double> day_summer;

    QVector<double> dt_transition;
    QVector<double> day_transition;
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void Add(sundial * a);
    void Paint();
    void Paint2();
private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

private:
    Ui::MainWindow *ui;
    sundial *arg;
};

#endif // MAINWINDOW_H
