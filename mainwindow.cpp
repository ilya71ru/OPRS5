#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::Add(sundial *a)
{
    arg=std::move(a);
}

void MainWindow::Paint()
{
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(arg->erg[2],  arg->erg[0]);
    ui->widget->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->widget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle,1));
    ui->widget->xAxis->setLabel("Z");
    ui->widget->yAxis->setLabel("X");
   ui->widget->xAxis->setRange(-20,20);
    ui->widget->yAxis->setRange(-20,20);


    ui->widget->replot();
}


void MainWindow::Paint2()
{
    ui->widget_2->clearGraphs();
    ui->widget_2->addGraph();
     ui->widget_2->addGraph();
      ui->widget_2->addGraph();
       ui->widget_2->addGraph();
    ui->widget_2->graph(0)->setData(day,  dt);
    if (ui->checkBox_2->isChecked()==true)
    {

    ui->widget_2->graph(1)->setData(day,  dt_between);
    ui->widget_2->graph(1)->setPen(QPen(QColor(Qt::red)));
    }
    if (ui->checkBox->isChecked()==true)
    {

        ui->widget_2->graph(2)->setData(day_summer,  dt_summer);
        ui->widget_2->graph(2)->setPen(QPen(QColor(Qt::green)));
    }
    if (ui->checkBox_3->isChecked()==true)
    {

        ui->widget_2->graph(3)->setData(day_transition,  dt_transition);
        ui->widget_2->graph(3)->setPen(QPen(QColor(Qt::yellow)));
    }
    ui->widget_2->xAxis->setLabel("day");
    ui->widget_2->yAxis->setLabel("time");
    ui->widget_2->xAxis->setRange(1,365);
    ui->widget_2->yAxis->setRange(0,24);
    ui->widget_2->replot();
}

void MainWindow::on_pushButton_clicked()
{
    ui->widget->xAxis->setRange((ui->lineEdit->text()).toInt(),ui->lineEdit_2->text().toInt());
    ui->widget->yAxis->setRange((ui->lineEdit_3->text()).toInt(),ui->lineEdit_4->text().toInt());
    ui->widget->replot();
}

void MainWindow::on_pushButton_2_clicked()
{

    TVector X0(6);
    X0.SetItem(0, -2.566123740124270e+7L);
    X0.SetItem(1, 1.339350231544666e+8L);
    X0.SetItem(2, 5.805149372446711e+7L);
    X0.SetItem(3, -29.83549561177192L);
    X0.SetItem(4, -4.846747552523134L);
    X0.SetItem(5, -2.100585886567924L);

    std::vector<double> month;
    month.push_back(31); //J
     month.push_back(28); //F
      month.push_back(31); //M
       month.push_back(30); //A
        month.push_back(31); //M
         month.push_back(30); //J
          month.push_back(31); //July
           month.push_back(31); //A
            month.push_back(30); //S
             month.push_back(31); //O
              month.push_back(30); //N
               month.push_back(31); //D
    int DYear=ui->lineEdit_8->text().toInt()-2019;
    int Month=ui->lineEdit_9->text().toInt();
    int Day=0;

    for (int i=0;i<Month-1;i++)
    {
        Day+=month.at(i);
    }

    Day=Day+ui->lineEdit_10->text().toInt();

    long double All=DYear*365*86400+Day*86400;

    sundial *Arg1=new sundial(0,All,All,X0);
    Arg1->SetGnomon(ui->lineEdit_5->text().toInt(),ui->lineEdit_6->text().toInt(),ui->lineEdit_7->text().toInt());
    TIntegrator * integrator=new TDormandPrince();
    integrator->setPrecision(1e-16);
    integrator->Run(Arg1);
    Arg1->SetTime(2019,1,1,0,0,0);
    TMatrix Result = Arg1->getResult();

    std::ofstream out;
    out.open("first.txt");
    if (out.is_open())
    {
           for (int i=0;i<Result.GetRowCount();i++)
                   {
                           for(int j=1;j<Result.GetColCount()-3;j++)
                           {
                               out<< std::to_string(Result(i,j));
                               out<<(" ");
                           }
                           out <<std::endl;
                   }
                   out.close();
    }

TVector X(6);
X[0]=Result(Result.GetRowHigh(),1);
X[1]=Result(Result.GetRowHigh(),2);
X[2]=Result(Result.GetRowHigh(),3);
X[3]=Result(Result.GetRowHigh(),4);
X[4]=Result(Result.GetRowHigh(),5);
X[5]=Result(Result.GetRowHigh(),6);

    sundial *Arg=new sundial(0, 86400, 60, X);
     Arg->SetGnomon(ui->lineEdit_5->text().toDouble()+ui->lineEdit_11->text().toDouble()/60.,ui->lineEdit_6->text().toDouble()+ui->lineEdit_12->text().toDouble()/60.,ui->lineEdit_7->text().toDouble());
     Arg->SetTime(ui->lineEdit_8->text().toInt(),ui->lineEdit_9->text().toInt(),ui->lineEdit_10->text().toInt(),0,0,0);

    integrator->Run(Arg);

     Arg->toRotate(0,86400);
     out.open("eee.txt");
     int index=0;
     if (out.is_open())
     {
                     for(auto i=Arg->erg[0].begin();i<Arg->erg[0].end();i++)
                                          {

                                            out<<(" ");
                                            out<<std::to_string(Arg->erg[0].at(index));
                                            out<<(" ");
                                            out<<std::to_string(Arg->erg[2].at(index));

                                            out <<std::endl;

                                            index++;
                                         }
                    out.close();
     }
     Add(Arg);
     Paint();
     delete integrator;
     delete Arg1;
     delete Arg;
}


\

void MainWindow::on_pushButton_3_clicked()
{
    int start=8;
    int end=20;

    TVector X(6);
    X.SetItem(0, -2.566123740124270e+7L);
    X.SetItem(1, 1.339350231544666e+8L);
    X.SetItem(2, 5.805149372446711e+7L);
    X.SetItem(3, -29.83549561177192L);
    X.SetItem(4, -4.846747552523134L);
    X.SetItem(5, -2.100585886567924L);
     //Clear data
     this->dt.clear();
     this->day.clear();
     this->dt_between.clear();
     this->dt_summer.clear();
     this->day_summer.clear();
     this->dt_transition.clear();
     this->day_transition.clear();
    int timezone=ui->comboBox_2->currentText().toInt();
    sundial *Arg2=new sundial(0,86400*365,60,X);
    Arg2->setTimeZone(timezone);
    Arg2->SetGnomon(ui->lineEdit_13->text().toDouble()+ui->lineEdit_15->text().toDouble()/60.,ui->lineEdit_14->text().toDouble()+ui->lineEdit_16->text().toDouble()/60.,ui->lineEdit_17->text().toDouble());
    Arg2->SetTime(2019,1,1,0,0,0);
    TIntegrator * integrator2=new TDormandPrince();
    integrator2->setPrecision(1e-16);
    integrator2->Run(Arg2);
    Arg2->toRotate(0,86400*365);
    for(auto i=0;i<365;i++)
    {
//between 8 and 20 hours
double time_2=0;
double time=Arg2->End[i]/3600.-Arg2->Start[i]/3600.;

if ((ui->checkBox_2->isChecked()==true))
{

    double max1=0;
    double max2=0;
    if ((Arg2->End[i])>end*3600) max2=end*3600; else max2=Arg2->End[i];
    if ((Arg2->Start[i])>start*3600) max1=Arg2->Start[i]; else max1=start*3600;
    time_2=(max2-max1)/3600.;

}

dt.push_back(time);
day.push_back(i+1);
dt_between.push_back(time_2);

}
    delete integrator2;
    delete Arg2;

    std::ofstream out;
    out.open("polden.txt");
    if (out.is_open())
    {
           for (int i=0;i<365;i++)
                   {

                           out<< std::to_string(i+1);
                           out<<(" ");
                           out<< std::to_string(Arg2->Polden[i]);
                           out <<std::endl;
                   }
                   out.close();
    }

    X.SetItem(0, -2.566123740124270e+7L);
    X.SetItem(1, 1.339350231544666e+8L);
    X.SetItem(2, 5.805149372446711e+7L);
    X.SetItem(3, -29.83549561177192L);
    X.SetItem(4, -4.846747552523134L);
    X.SetItem(5, -2.100585886567924L);
if (ui->checkBox->isChecked()==true)
{
    int timezone_summer=ui->comboBox_2->currentText().toInt()+1;
    sundial *Arg3=new sundial(0,86400*365,60,X);
    Arg3->setTimeZone(timezone_summer);
    Arg3->SetGnomon(ui->lineEdit_13->text().toDouble()+ui->lineEdit_15->text().toDouble()/60.,ui->lineEdit_14->text().toDouble()+ui->lineEdit_16->text().toDouble()/60.,ui->lineEdit_17->text().toDouble());
    Arg3->SetTime(2019,1,1,0,0,0);
    TIntegrator * integrator3=new TDormandPrince();
    integrator3->setPrecision(1e-16);
    integrator3->Run(Arg3);
    Arg3->toRotate(0,86400*365);
    double oldtime=0;
    for(auto i=0;i<365;i++)
    {
//between 8 and 20 hours
double time_2=0;

double max1=0;
double max2=0;
if ((Arg2->End[i])>end*3600) max2=end*3600; else max2=Arg2->End[i];
if ((Arg2->Start[i])>start*3600) max1=Arg2->Start[i]; else max1=start*3600;
time_2=(max2-max1)/3600.;


dt_summer.push_back(time_2);
day_summer.push_back(i+1);

    }

    std::ofstream out;
    out.open("PoldenSummer.txt");
    if (out.is_open())
    {
           for (int i=0;i<365;i++)
                   {

                           out<< std::to_string(i+1);
                           out<<(" ");
                           out<< std::to_string(Arg2->Polden[i]);
                           out <<std::endl;
                   }
                   out.close();
    }
delete integrator3;
delete Arg3;

}

X.SetItem(0, -2.566123740124270e+7L);
X.SetItem(1, 1.339350231544666e+8L);
X.SetItem(2, 5.805149372446711e+7L);
X.SetItem(3, -29.83549561177192L);
X.SetItem(4, -4.846747552523134L);
X.SetItem(5, -2.100585886567924L);
  timezone=ui->comboBox_2->currentText().toInt();
if (ui->checkBox_3->isChecked()==true)
{
    sundial *Arg4=new sundial(0,86400*365,60,X);
int timezone_summer=ui->comboBox_2->currentText().toInt()+1;

    Arg4->SetGnomon(ui->lineEdit_13->text().toDouble()+ui->lineEdit_15->text().toDouble()/60.,ui->lineEdit_14->text().toDouble()+ui->lineEdit_16->text().toDouble()/60.,ui->lineEdit_17->text().toDouble());
    Arg4->SetTime(2019,1,1,0,0,0);
    TIntegrator * integrator4=new TDormandPrince();
    integrator4->setPrecision(1e-16);
    integrator4->Run(Arg4);
    Arg4->setTimeZone(timezone);
    Arg4->toRotate(0,86400*89+2*3600);
    Arg4->setTimeZone(timezone+1);
    Arg4->toRotate(86400*89+2*3600,86400*301+3*3600);
    Arg4->setTimeZone(timezone);
    Arg4->toRotate(86400*301+3*3600,86400*365);

    for(auto i=0;i<365;i++)
    {
//between 8 and 20 hours
double time_2=0;

double max1=0;
double max2=0;
if ((Arg2->End[i])>end*3600) max2=end*3600; else max2=Arg2->End[i];
if ((Arg2->Start[i])>start*3600) max1=Arg2->Start[i]; else max1=start*3600;
time_2=(max2-max1)/3600.;

dt_transition.push_back(time_2);
day_transition.push_back(i);
}

    std::ofstream out;
    out.open("PoldenBeetween.txt");
    if (out.is_open())
    {
           for (int i=0;i<365;i++)
                   {

                           out<< std::to_string(i+1);
                           out<<(" ");
                           out<< std::to_string(Arg2->Polden[i]);
                           out <<std::endl;
                   }
                   out.close();
    }

    delete integrator4;
    delete Arg4;

}
FILE *gp=popen("gnuplot -persist","w");
if (gp==nullptr)
{
    printf("error opening pipe to GNU");
    exit(1);
}
fprintf(gp, "plot \"Polden.txt\" w l,\"PoldenSummer.txt\" w l,\"PoldenBeetween.txt\" w l \n");
pclose(gp);
Paint2();
}
