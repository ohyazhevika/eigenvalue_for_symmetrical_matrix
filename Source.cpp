#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include "Matrix.h"
using namespace std;

double Rand(int range)
{
	double ans = (double)(rand() % (2 * abs(range) + 1) - abs(range));
	ans += 1.0*(rand() % 1000) / 1000;
	return ans;
}

double EuclideanNorm(Matrix M, int n) //m - вектор-столбцом должен быть
{
	double res = 0;
	double comp;
	for (int i = 1; i <= n; i++)
	{
		comp=M(i, 1);
		res += comp * comp;
	}
	return sqrt(res);
}

double FirstNorm(Matrix M, int n) //перва€ норма дл€ вектора (не дл€ матрицы)
{
	double res = 0;
	double fcomp;
	for (int i = 1; i <= n; i++)
	{
		fcomp = fabs(M(i, 1));
		if (fcomp > res)
			res = fcomp;
	}
	return res;
}
int main()
{                                                                               
	srand(time(NULL));
	setlocale(LC_ALL, "Russian");
	int N = 50; //size of matrix
	int range = 50;
	double epsilon = 0.00000001;
	double 
		avg_r = 0,				// —р. мера точности r
		avg_lambda_rel_err = 0, //—р. оценка точности соб. значений
		avg_vector_rel_err = 0; //—р. оценка точности соб. векторов
	int avg_iter_cnt = 0;		//—р. число итераций
	int numberOfTests = 10;
	for (int i=0; i<numberOfTests; i++)
	{
		double
			* eigenvalues = new double[N],  //соб значени€ матрицы
			* w = new double[N],
			* e = new double[N],	// единичный вектор
			 wlength = 0;
		int maxLamInd = 0;			//индекс максимального по модулю собств значени€ в eigenvalues
		int maxLam2Ind = 0;			//второй максимальный по модулю
		for (int i = 0; i < N; i++)
		{
			eigenvalues[i] = Rand(range);
			if (fabs(eigenvalues[i]) > fabs(eigenvalues[maxLamInd]))
			{
				maxLam2Ind = maxLamInd;
				maxLamInd = i;
			}
			else
				if (fabs(eigenvalues[i]) > fabs(eigenvalues[maxLam2Ind]))
					maxLam2Ind = i;
			w[i] = Rand(range);
			e[i] = 1;
			wlength += w[i] * w[i];
		}
		wlength = 1 / sqrt(wlength);
		for (int i = 0; i < N; i++)
			w[i] *= wlength;
		Matrix Lambda(N, eigenvalues);	//диагональна€ матрица с соб значени€ми на диагнонали
		Matrix W(N, 1, w);	//вектор w
		Matrix Wt(1, N, w);	//транспонированный вектор
		Matrix E(N, e);		//единична€ матрица
		Matrix H = E - (W * Wt) * 2;
		Matrix A = H * Lambda * H;
		double* xlist = new double[N];
		double* xlist2 = new double[N];
		for (int i = 0; i < N; i++)
		{
			xlist[i] = H(i + 1, maxLamInd + 1); //компоненты реального собственного вектора, отвечающего макс. соб значению
			xlist2[i] = H(i + 1, maxLam2Ind + 1);
		}
		Matrix XLn(N, 1, xlist);	//соб вектор (столбец), отвечающий макс. соб. значению
		Matrix XLnT(1, N, xlist);	//соб вектор (строка), отвечающий макс. соб. значению (транспонированный столбец)
		Matrix XLn2(N, 1, xlist2);	//соб вектор (столбец), отвечающий 2-му макс. соб значению
		double lambdaN = eigenvalues[maxLamInd];	//максимальное по модулю собственное значение
		double lambdaN2 = eigenvalues[maxLam2Ind];	//второе макс. по модулю собственное значение
		Matrix A1 = A - (XLn * XLnT) * lambdaN;
		Matrix X(N, 1, e);	//первое приближение X(0) (столбец)
		Matrix V = X * (1 / EuclideanNorm(X, N));
		X = A1 * V;
		double Sigma = 0;
		for (int i = 1; i <= N; i++)
			Sigma += V(i, 1) * X(i, 1);
		double prev_sigma = Sigma;
		Matrix prev_X = X;
		int numberOfIter = 0;
		do
		{
			numberOfIter++;
			prev_sigma = Sigma;
			prev_X = X;
			V = X * (1 / EuclideanNorm(X, N));
			X = A1 * V;
			Sigma = 0;
			for (int i = 1; i <= N; i++)
				Sigma += V(i, 1) * X(i, 1);
		} while (fabs(prev_sigma - Sigma) > epsilon && (numberOfIter < 10000));
		double lambda_rel_error = (lambdaN2 == 0) ? fabs(Sigma) : fabs((Sigma - lambdaN2) / lambdaN2);
		double vector_rel_error = 0;
		for (int i = 1; i <= N; i++)
		{
			double
				calc = V(i, 1),
				real = XLn2(i, 1),
				err = (real == 0) ? fabs(calc) : fabs((real - calc) / real);
			vector_rel_error = max(vector_rel_error, err);

		}
		avg_iter_cnt += numberOfIter;
		avg_lambda_rel_err += lambda_rel_error;
		avg_r += FirstNorm(A * V - V * Sigma, N);
		avg_vector_rel_err += vector_rel_error;
	}
	avg_iter_cnt /= numberOfTests;
	avg_lambda_rel_err /= numberOfTests;
	avg_r /= numberOfTests;
	avg_vector_rel_err /= numberOfTests;
	cout
		<< "—редн€ мера точности r = " << avg_r << endl
		<< "—реднее число итераций = " << avg_iter_cnt << endl
		<< "—редн€€ оценка точности соб. значений = " << avg_lambda_rel_err << endl
		<< "—редн€€ оценка точности соб. векторов = " << avg_vector_rel_err << endl;
	system("pause");
}