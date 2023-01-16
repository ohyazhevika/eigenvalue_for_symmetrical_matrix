#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>

using namespace std;


class Matrix
{
protected:
	int rows; //Количество строк матрицы
	int cols; //Количество столбцов матрицы
	double** cells;
	void AllocateCells(int, int); // Выделение памяти под матрицу
	void FreeCells();
public:
	Matrix() : rows(0), cols(0), cells(nullptr) {}	// Конструктор по умолчанию
	Matrix(const Matrix&);					// Конструктор копирования
	Matrix(int, int);							// Конструктор нулевой матрицы
	Matrix(int, double*);						// Конструктор квадратнрой диагональной матрицы из списка
	Matrix(int, int, double*);				//Конструктор матрицы из списка (целиком!)
	~Matrix();								// Деструктор

	double& operator()(int i, int j) { return cells[i - 1][j - 1]; }

	Matrix& operator = (const Matrix&);		// Перегрузка оператора присваивания
	Matrix  operator + (const Matrix&);		// Сложение матриц
	Matrix  operator - (const Matrix&);		// Вычитание матриц
	Matrix  operator * (const Matrix&);		// Умножение матриц
	Matrix operator *(double);

	Matrix Transpone();						//Транспонирование матрицы

	friend istream& operator >> (istream&, Matrix&);			// Перегрузка оператора >> для ввода матрицы
	friend ostream& operator << (ostream&, Matrix&);	// Перегрузка оператора << для вывода матрицы
	int GetRows() { return rows; }
	int GetCols() { return cols; }
};



Matrix::Matrix(const Matrix& M)
{
	AllocateCells(M.rows, M.cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			cells[i][j] = M.cells[i][j];
}

Matrix::Matrix(int Rows, int Cols)
{
	AllocateCells(Rows, Cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			cells[i][j] = 0;
}


Matrix::Matrix (int Rows, double* list)
{
	AllocateCells(Rows, Rows);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < rows; j++)
		{
			if (i != j)
				cells[i][j] = 0;
			else
				cells[i][i] = list[i];
		}	
}

Matrix::Matrix(int Rows, int Cols, double* list)
{
	AllocateCells(Rows, Cols);
	int count = 0;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; count++, j++)
			cells[i][j] = list[count];
}

Matrix::~Matrix()
{
	FreeCells();
}


Matrix& Matrix::operator=(const Matrix& M)
{
	if (!(rows == M.rows && cols == M.cols))
	{
		FreeCells();
		AllocateCells(M.rows, M.cols);
	}
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			cells[i][j] = M.cells[i][j];
	return *this;
}

Matrix Matrix::operator+(const Matrix& M)
{
	Matrix res(*this);
	if (rows == M.rows && cols == M.cols)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				res.cells[i][j] += M.cells[i][j];
	}
	return res;
}


Matrix Matrix::operator-(const Matrix& M)
{
	Matrix res(*this);
	if (rows == M.rows && cols == M.cols)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				res.cells[i][j] -= M.cells[i][j];
	}
	return res;
}


Matrix Matrix::operator*(const Matrix& M)
{
	Matrix res(rows, M.cols);
	if (cols == M.rows)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < M.cols; j++)
			{
				double summ = 0;
				for (int k = 0; k < cols; k++)
					summ += cells[i][k] * M.cells[k][j];
				res.cells[i][j] = summ;
			}
		}
	}
	return res;
}

Matrix Matrix::operator*(double a)
{
	Matrix res = *this;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			res.cells[i][j] *= a;
	return res;
}


void Matrix::AllocateCells(int Rows, int Cols)
{
	cells = new double * [Rows];
	for (int i = 0; i < Rows; i++)
		cells[i] = new double[Cols];
	rows = Rows;
	cols = Cols;
}


void Matrix::FreeCells()
{
	for (int i = 0; i < rows; i++)
		delete cells[i];
	delete cells;
	rows = cols = 0;
}

istream& operator >> (istream& fi, Matrix& M)
{
	for (int i = 1; i <= M.GetRows(); i++)
		for (int j = 1; j <= M.GetCols(); j++)
			fi >> M(i,j);
	return fi;
}

ostream& operator << (ostream& fo, Matrix& M)
{
	for (int i = 1; i <= M.GetRows(); i++)
	{
		//fo << "  ";
		for (int j = 1; j <= M.GetCols(); j++)
			fo << M(i, j)
			   << " \t";
		fo << endl;
	}
	return fo;
}

#endif MATRIX_H
