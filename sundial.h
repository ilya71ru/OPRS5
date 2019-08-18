#ifndef SUNDIAL_H
#define SUNDIAL_H
#include "model.h"
#include "spaceformule.h"
#include <QVector>
class sundial : public TModel
{
protected:
    //TVector X;

    const double mu_s = 132712.43994e+6;
    const double omega= 7.292115e-5;

    const long double R_=6371.3L;
public:
    TVector RE_;
    TVector RSH;
    TVector R0;
    TVector Rg;
    TVector Re0;
    TVector R;
    long double S;
    long double S_G0;
    long double height;
    long double L;
    long double Fi;

    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;

    int index=0;
    double timezone=0;
    bool daytime = false; //is it day or night
    QVector <double> erg[5];
    long double Polden[365];

    long double Start[365];
    long double End[365];
    sundial(long double t0, long double t1, long double SamplingIncrement,const TVector& X);
    TVector getRight( const TVector& X, long double t);
   void SetTime(int year, int month, int day, int hour, int minute, int second);
   void SetGnomon(long double L,long double Fi,long double H);
   void toRotate(long double t0, long double t1);
   void setTimeZone(double t);
};

#endif // SUNDIAL_H
