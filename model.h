//---------------------------------------------------------------------------

#ifndef modelH
#define modelH

// Подключение заголовочного файла библиотеки векторных и матричных операций
#include "MyLinearAlgebra.h"
#include <math.h>

//---------------------------------------------------------------------------
// Базовый класс модели для интегратора
using namespace MyLinearAlgebra;

class TModel
{
    protected:
        // Начальные условия
	    TVector X0;
        // Требуемый от интегратора интервал выдачи результатов
        long double SamplingIncrement;
        // Начало и окончание времени интегрирования
        long double t0, t1;
        // Хранилище результатов
        TMatrix Result;
        // Счётчик строк в матрице результатов
		int N;
		
    public:
        // Общий конструктор моделей - инициализация параметров по умолчанию
        TModel()
            : SamplingIncrement( 1e-1 )
            , t0( 0 )
            , t1( 1 )
			, N( 0 )
        {}
        // Абстрактная перегружаемая функция правых частей ДУ (X - вектор состояния, t - независимый аргумент)
        virtual TVector getRight( const TVector& X, long double t ) = 0;
        // Получение начальных условий
        inline TVector getInitialConditions() const { return X0; }
        // Порядок системы - по размерности вектора начальных условий
        inline int getOrder() const { return X0.GetSize(); }

        // Интервал выдачи результатов
        inline long double getSamplingIncrement() const { return SamplingIncrement; }

        // Управление временным интервалом интегрирования
        inline long double getT0() const { return t0; }
        inline long double getT1() const { return t1; }
       
		// Получение матрицы результатов
		inline TMatrix getResult() const { return Result; }
        // Запись результатов (в этом методе в наследнике допустимо организовать запись в файл 
		// вместо или вместе с наполнением матрицы результатов)
        virtual void addResult( const TVector& X, long double &t );
		// Очистка матрицы результатов (или файла с результатами)
		virtual void clearResult();
		// Подготовка матрицы результатов для более эффективного выделения памяти
		virtual void prepareResult();
		virtual ~TModel(){};
};

//---------------------------------------------------------------------------
class TIntegrator
{
    protected:
        // Максимальная погрешность на шаге
        long double Eps;
    public:
		// Базовый конструктор
        TIntegrator() : Eps( 1e-8 ) {}
        inline void setPrecision( long double Eps ) { this->Eps = Eps; }
        inline long double getPrecision() const { return Eps; }
		// Абстрактный метод, реализующий процедуру численного интегрирования и возвращающий глобальную погрешность вычислений
        virtual void Run(TModel* Model)=0;
        virtual ~TIntegrator(){}
};

//---------------------------------------------------------------------------

class TDormandPrince : public TIntegrator
{
    private:
        // Коэффициенты a,b,c
    const long double c[7] = { 0, 1./5, 3./10, 4./5, 8./9, 1., 1. };
	const long double a[7][6] = {
	    { 0. },
	    { 1./5 },
	    { 3./40, 9./40 },
	    { 44./45, -56./15, 32./9 },
	    { 19372./6561, -25360./2187, 64448./6561, -212./729 },
	    { 9017./3168, -355./33, 46732./5247, 49./176, -5103./18656 },
	    { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84 }
	};
	const long double b1[7] = { 35./384, 0., 500./1113, 125./192, -2187./6784, 11./84, 0 };
	const long double b2[7] = { 5179./57600, 0., 7571./16695, 393./640, -92097./339200, 187./2100, 1./40 };

        // Коэффициенты K[i]
        TVector K[7];
		// Машинный нуль
		long double u;
    public:
        TDormandPrince();

        virtual void Run(TModel* Model);

};


#endif
 
