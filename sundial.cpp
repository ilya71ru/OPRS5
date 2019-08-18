#include "sundial.h"
#include <cmath>
using namespace MyLinearAlgebra;


sundial::sundial(long double t0, long double t1, long double SamplingIncrement,const TVector& X)
    : TModel()
{
    this->t0=t0;
    this->t1=t1;
    this->SamplingIncrement=SamplingIncrement;
    RE_=TVector(3);
    RSH=TVector(3);//то что нужно найти
    R0=TVector(3);
    Rg=TVector(3);
    R=TVector(3);
    Re0=TVector(3);
    for(int i=0;i<365;i++)
      {
        Start[i]=-1;
        End[i]=-1;
        Polden[i]=0;
      }
    X0=X;
}

TVector sundial::getRight( const TVector& X, long double t) {
    TVector Y(6);
    for (int i=0; i!=3; ++i)
        Y.SetItem(i, X[i+3]);
    TVector X_Coordin = X.Concat(0,2);

    TVector X_Dif = -mu_s*X_Coordin*(1/(pow(X_Coordin.length(),3)));
    for (int i=3; i!=6; ++i)
        Y.SetItem(i, X_Dif[i-3]);  



    return Y;
}



void sundial::SetTime(int year, int month, int day, int hour, int minute, int second)
{
    this->year=year;
    this->month=month;
    this->day=day;
    this->hour=hour;
    this->minute=minute;
    this->second=second;
    S_G0=2*M_PI/86400L*fmod(FSpace::SiderialTime(year, month,day,hour, minute, second),86400L);
}

void sundial::SetGnomon(long double L, long double Fi, long double H)
{
     this->L=L*M_PI/180;
     this->Fi=Fi*M_PI/180;
     this->height=H/1000L;
}

void sundial::toRotate(long double t0, long double t1)
{
    long double P=100000;
    int key=index;
    TVector X_Coordin(3);
    while(t0<t1)
    {

        for(int i=0;i<3;i++)
            X_Coordin[i]=Result(t0/SamplingIncrement,i+1);


        S=S_G0+omega*t0+L;


        for(int i = 0; i < 3; i++) Re0[i] = X_Coordin[i]/X_Coordin.length()/*ro*/;

         R[0] = R_*cos(Fi)*cos(S);
         R[1] = R_*cos(Fi)*sin(S);
         R[2] = R_*sin(Fi);

         R0=R*(1/R.length());

        Rg =height*R0;

        RE_=-height/(Re0*R0)*Re0;

        RSH=RE_+Rg;

        TMatrix A=FSpace::intoA(Fi,S);

        RSH=A*RSH;

        long double alfa=X_Coordin*Rg;



        if (alfa < 0) {daytime = true;
            if (Start[index]==-1) Start[index]=t0+timezone*3600-index*86400;
          }
            else {
            if ((Start[index]!=-1)  && (End[index]==-1)) {End[index]=t0+timezone*3600-index*86400;  index++;}
            daytime=false;

         }

         double c=1000*sqrt(RSH[0]*RSH[0]+RSH[2]*RSH[2]);
        if ((P>c)&&(key==index)&&(daytime==true))
          {
            P=c;
          }
        else if ((daytime==true)&&(key==index)){
            if((t0+timezone*3600-(int)(t0/86400)*86400)/3600<=24)
              Polden[index]=(t0+timezone*3600-(int)(t0/86400)*86400)/3600;
            else Polden[index]=Polden[index-1];
              key++;
          }
        else if (daytime==false)
          {
            P=100000;
          }


        if (daytime == true){
        //    Time[(int)(t0/86400)].push_back(t0-(int)(t0/86400)*86400);



            erg[3].push_back(t0);

            double c=1000*sqrt(RSH[0]*RSH[0]+RSH[2]*RSH[2]);
            if (c<20*height*1000)
            {     
         erg[0].push_back(RSH[0]*1000);
          erg[1].push_back(RSH[1]*1000);
           erg[2].push_back(RSH[2]*1000);
             erg[4].push_back(c);
            }
         }
         daytime=false;

         t0+=SamplingIncrement;
    }
}

void sundial::setTimeZone(double t)
{
  timezone=t;
}
