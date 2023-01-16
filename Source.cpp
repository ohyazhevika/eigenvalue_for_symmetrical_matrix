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

double EuclideanNorm(Matrix M, int n) //m - ������-�������� ������ ����
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

double FirstNorm(Matrix M, int n) //������ ����� ��� ������� (�� ��� �������)
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
		avg_r = 0,				// ��. ���� �������� r
		avg_lambda_rel_err = 0, //��. ������ �������� ���. ��������
		avg_vector_rel_err = 0; //��. ������ �������� ���. ��������
	int avg_iter_cnt = 0;		//��. ����� ��������
	int numberOfTests = 10;
	for (int i=0; i<numberOfTests; i++)
	{
		double
			* eigenvalues = new double[N],  //��� �������� �������
			* w = new double[N],
			* e = new double[N],	// ��������� ������
			 wlength = 0;
		int maxLamInd = 0;			//������ ������������� �� ������ ������ �������� � eigenvalues
		int maxLam2Ind = 0;			//������ ������������ �� ������
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
		Matrix Lambda(N, eigenvalues);	//������������ ������� � ��� ���������� �� ����������
		Matrix W(N, 1, w);	//������ w
		Matrix Wt(1, N, w);	//����������������� ������
		Matrix E(N, e);		//��������� �������
		Matrix H = E - (W * Wt) * 2;
		Matrix A = H * Lambda * H;
		double* xlist = new double[N];
		double* xlist2 = new double[N];
		for (int i = 0; i < N; i++)
		{
			xlist[i] = H(i + 1, maxLamInd + 1); //���������� ��������� ������������ �������, ����������� ����. ��� ��������
			xlist2[i] = H(i + 1, maxLam2Ind + 1);
		}
		Matrix XLn(N, 1, xlist);	//��� ������ (�������), ���������� ����. ���. ��������
		Matrix XLnT(1, N, xlist);	//��� ������ (������), ���������� ����. ���. �������� (����������������� �������)
		Matrix XLn2(N, 1, xlist2);	//��� ������ (�������), ���������� 2-�� ����. ��� ��������
		double lambdaN = eigenvalues[maxLamInd];	//������������ �� ������ ����������� ��������
		double lambdaN2 = eigenvalues[maxLam2Ind];	//������ ����. �� ������ ����������� ��������
		Matrix A1 = A - (XLn * XLnT) * lambdaN;
		Matrix X(N, 1, e);	//������ ����������� X(0) (�������)
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
		<< "������ ���� �������� r = " << avg_r << endl
		<< "������� ����� �������� = " << avg_iter_cnt << endl
		<< "������� ������ �������� ���. �������� = " << avg_lambda_rel_err << endl
		<< "������� ������ �������� ���. �������� = " << avg_vector_rel_err << endl;
	system("pause");
}