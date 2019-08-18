#include <QApplication>
#include <iostream>
#include <cmath>
#include "spaceformule.h"
#include "sundial.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "mainwindow.h"
//using namespace std;



int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow window;

    window.show();




    return a.exec();
}


