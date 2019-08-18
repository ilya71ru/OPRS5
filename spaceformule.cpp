#include "spaceformule.h"

namespace FSpace {

double JulianDate(int year, int month, int day, int hour, int minute, int second) {
    int a = (14 - month) / 12;
    int Y = year + 4800 - a;
    int M = month + 12 * a - 3;
    int JDN = day + (153 * M + 2) / 5 + 365 * Y+ Y / 4 - Y / 100 + Y / 400 -32045;
    double JD = JDN + (hour - 12)/24.0 + minute/1440.0 + second/86400.0;
    return JD;
}

double SiderialTime(int year, int month, int day, int hour, int minute, int second){
    int d = (int)(JulianDate(year, month, day, hour, minute, second) - 2451544.5);
    double t = d / 36525.0;
    double S_G_0 = 24110.584841 + 8640184.812866 * t + 0.093104 * pow(t, 2.) - 6.2e-6 * pow(t, 3.);
    return S_G_0;
}

TVector intoGSK(TVector R)
{
     long double RE=6371.3L;
     TVector r(3);
     r[0] = RE*cos(R[1])*cos(R[0]);
     r[1] = RE*cos(R[1])*sin(R[0]);
     r[2] = RE*sin(R[1]);
     return r;
}

TMatrix intoA(double Fi,double S)
{
    TMatrix A(3,3);
    A(0,0)=-sin(Fi)*cos(S);
    A(0,1)=-sin(Fi)*sin(S);
    A(0,2)=cos(Fi);
    A(1,0)=cos(Fi)*cos(S);
    A(1,1)=cos(Fi)*sin(S);
    A(1,2)=-sin(Fi);//?Fi
    A(2,0)=-sin(S);
    A(2,1)=cos(S);
    A(2,2)=0;
    return A;
}
}
