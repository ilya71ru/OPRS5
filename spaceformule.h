#ifndef SPACEFORMULE_H
#define SPACEFORMULE_H
#include <cmath>
#include "MyLinearAlgebra.h"
using namespace MyLinearAlgebra;
namespace FSpace{

double JulianDate(int year, int month, int day, int hour, int minute, int second);
double SiderialTime(int year, int month, int day, int hour, int minute, int second);
TVector intoGSK(TVector R);
TMatrix intoA(double Fi,double S);

}
#endif // SPACEFORMULE_H
