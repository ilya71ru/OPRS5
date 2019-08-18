/*
 * TVector.h
 *
 *  Created on: 31 дек. 2018 г.
 *      Author: Maks
 */

#ifndef MYLINEARALGEBRA_H_
#define MYLINEARALGEBRA_H_
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string.h>

namespace MyLinearAlgebra {
class TMatrix;
class TQuaternion;
class TVector {
protected:
int n;
long double *data;
public:
	TVector() noexcept;
	TVector(int n) noexcept;

	TVector(long double a[],int n);

	TVector(const TVector& value);
	virtual ~TVector();
	inline long double GetItem(int i)  {return data[i];}
	inline void SetItem(int i,long double value)  {data[i]=value;}
    inline int GetSize() const noexcept  {return n;}
	inline int GetHight() const noexcept {return n-1;}
	inline long double& operator[](int i){return data[i];}
	inline const long double& operator[](int i) const {return data[i];}
	TVector& operator=(const TVector& value);
	TVector Clone();
	void Resize(int n);
	TVector operator-()const;
	TVector operator-(const TVector& a) const;
	TVector operator+(const TVector& a) const;
	TVector operator *(long double a) const;
	long double operator *(const TVector& a) const;
	TVector operator *(const TMatrix& a) const;
	TQuaternion operator *(const TQuaternion& a) const;
	TVector operator ^(const TVector& a) const;// векторное умножение
	friend TVector operator *(long double value,const TVector& a);
	long double length() const noexcept ; // получение модуля вектора
	TVector norm() noexcept ; //нормировка вектора
	void Print();
	TVector rotateByRodrigFormule(TVector& e, long double phi) const;
	TVector rotateByQuaternion(TQuaternion& Q) const;
	TVector toreturnKrilov(TQuaternion& Q) const;
	TMatrix inMatrix() const;// перевод в матрицу
    TVector Concat(int i, int j) const;// разрез вектора от i до j
	TVector Clip(TVector a); //склеивание векторов

};

class TMatrix
{
protected:

// Размерность матрицы (число строк и столбцов)

int n, m;

// Элементы матрицы

long double **data;

public:

// Конструктор по умолчанию

TMatrix() noexcept ;

// Конструктор с заданной размерностью

TMatrix(int n, int m) noexcept ;

// Конструктор копий

TMatrix(const TMatrix& rvalue);

TMatrix(long double ** A, int n, int m);

void Print() ;

// Оператор присваивания

TMatrix& operator = (const TMatrix& rvalue);

// Деструктор

virtual ~TMatrix();

// Функция получения количества строк

inline int GetRowCount() const noexcept { return n; }

// Функция получения кол-ва столбцов

inline int GetColCount() const noexcept  { return m; }

// Функция получения индекса последней строки

inline int GetRowHigh() const noexcept { return n-1; }

// Функция получения индекса последнего столбца

inline int GetColHigh() const noexcept { return m-1; }

// Функция задания размерности

void Resize(int n, int m);

// Оператор доступа к элементам матрицы

inline long double& operator()(int i, int j) { return data[i][j]; }

// Оператор константного доступа к элементам вектора

inline const long double& operator()(int i, int j) const { return data[i][j]; }

// Оператор - унарный минус

TMatrix operator - () const;

// Оператор вычитания матриц

TMatrix operator - (const TMatrix& arg) const;

// Оператор сложения матриц

TMatrix operator + (const TMatrix& arg) const;

// Оператор умножения матрицы на число

TMatrix operator * (long double arg) const;

// Оператор умножения матриц

TMatrix operator * (const TMatrix& arg) const;

// Оператор умножения матрицы на вектор

TVector operator * (const TVector& arg) const;

// Дружественная функция - оператор умножения числа на матрицу

friend TMatrix operator * (long double lvalue, const TMatrix& rvalue);

// Оператор обращения матриц (метод Гаусса)

TMatrix operator ! () const;

// Функция вычисления детерминанта

long double det() const;

// Функция транспонирования

TMatrix t() const noexcept ;

// Функция формирования единичной матрицы

static TMatrix E(int n) noexcept ;

// Функция перестановки строк

TMatrix swapRows(int i, int j);
TMatrix swapCols(int i, int j);
void toText(std::string name);
TVector inVector() const;

};
class TQuaternion
{
	// Скалярная часть

	long double q0;

	// Векторная часть

	TVector Q;

	public:

	// Конструктор по умолчанию

	TQuaternion();

	// Конструктор по компонентам кватерниона

	TQuaternion(long double l0, long double l1, long double l2, long double l3);

	// Конструктор по углу поворота (рад.) и оси вращения

	TQuaternion(long double phi, const TVector& e);

	// Конструктор копирования

	TQuaternion(const TQuaternion& rvalue);

	// Производящая функция для создания кватерниона по углам Крылова

	TQuaternion(long double yaw, long double pitch,long double roll);

	// Оператор присваивания

	TQuaternion& operator = (const TQuaternion& rvalue);

	// Оператор вычитания кватернионов

	TQuaternion operator - (const TQuaternion& arg) const;

	// Оператор сложения кватернионов

	TQuaternion operator + (const TQuaternion& arg) const;

	// Оператор умножения кватернионов

	TQuaternion operator * (const TQuaternion& arg) const;

	// Оператор умножения кватерниона на вектор

	TQuaternion operator * (const TVector& arg) const;

	// Оператор умножения кватерниона на число

	TQuaternion operator * (long double arg) const;

	// Дружественная функция - оператор умножения числа на кватернион

	friend TQuaternion operator * (long double lvalue, const TQuaternion & rvalue);

	// Оператор обращения кватерниона

	TQuaternion operator ! () const;

	// Доступ к скалярной части

	inline long double scal() const noexcept { return q0; }

	// Доступ к векторной части

	inline TVector vect() const noexcept { return Q; }

	// Функция нормирования кватерниона

	TQuaternion norm() noexcept ;

	//функция получения модуля кватерниона

	long double Length() const noexcept ;

	// Функция получения сопряженного кватерниона

	TQuaternion conj() const;

	// Функция получения матрицы вращения из компонентов кватерниона

	TMatrix toRotateMatrix() ;

	//функция возврата матрицы в кватернион

	TQuaternion toReturnQ(const TMatrix& A) const;

    void Print();
};
} /* namespace MyLinearAlgebra */

#endif /* TVECTOR_H_ */
