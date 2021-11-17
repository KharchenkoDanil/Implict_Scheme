#define _CRT_SECURE_NO_WARNINGS 
#define _USE_MATH_DEFINES
#include <stdio.h> 
#include <conio.h> 
#include <iostream>
#include <fstream>
#include <locale.h> 
#include <math.h> 
#include <cmath>
using namespace std;

#define pi 3.141592653589793238
#define l pi 	//права€ граница по пространству
#define T 2. 	//права€ граница по времени
#define N 20	//число узлов по пространству
#define nu 2.
#define K 100	
#define sigma 0.25

double *x = new double[N + 1];

double *_U = new double[N + 1];
double *U = new double[N + 1];
double *U_ = new double[N + 1];

double *d = new double[N + 1];


void nodes(double x[], double n)  // задаЄм равномерное разбиение
{
	for (int i = 0; i <= n; i++)
	{
		x[i] = i * l / n;
	}
}
double time(double r) { //определ€ем количество узлов по времени
	double y; 
	double s;
	y = modf(T / r, &s);
	return s;
 }

double u(double x, double t) // определ€ем функцию u(x)
{
	double z;
	z = cos(x)*exp(2. * t);
	return z;
}

double f(double x, double t) // определ€ем функцию f(x)
{
	double z;
	z = 4 * cos(x)*exp(2 * t);
	return z;
}

double u0(double x) //начальное условие при t=0
{
	double z;
	z = cos(x);
	return z;
}

double yO(double t) //краевое условие на левой границе при x=0
{
	double z;
	z = exp(2 * t);
	return z;
}

double yl(double t) //краевое условие на правой границе при x=l
{
	double z;
	z = -exp(2 * t);
	return z;
}

void solution(double *U_, double x[], double t[], double d[], int k, int n) 
{

	double *a = new double[n + 1];
	double *b = new double[n + 1];
	double *c = new double[n + 1];
	double h = x[1] - x[0];
	double r = t[1] - t[0];
	c[0] = -1.;
	b[0] = 0;
	a[n] = 0;
	c[n] = -1.;


	for (int i = 1; i <= n - 1; i++)
	{
		a[i] = nu*sigma / pow(h, 2.);
		b[i] = nu*sigma / pow(h, 2.);
		c[i] = 1 / (2 * r) + 2 * nu*sigma / pow(h, 2.);
	}

	double *alpha = new double[n + 1];
	double *betta = new double[n + 1];

	alpha[0] = b[0] / c[0];
	betta[0] = -d[0] / c[0];

	for (int j = 1; j <= n - 1; j++)
	{
		alpha[j] = b[j] / (c[j] - alpha[j - 1] * a[j]);
		betta[j] = (betta[j - 1] * a[j] - d[j]) / (c[j] - alpha[j - 1] * a[j]);
	}

	U_[n] = (betta[n - 1] * a[n] - d[n]) / (c[n] - alpha[n - 1] * a[n]);

	for (int i = n; i >= 1; i--)
	{
		U_[i - 1] = alpha[i - 1] * U_[i] + betta[i - 1];
	}


}

double max(double U[], double x[], double t[], int k, double n) //задаЄм норму
{
	double s = fabs(U[0] - u(x[0], t[k]));
	for (int i = 1; i <= n; i++)
	{

		if (fabs(U[i] - u(x[i], t[k])) >= s)
			s = fabs(U[i] - u(x[i], t[k]));
	}
	return s;
}

double maxEps(double *eps, int m)
{
	double s = eps[0];
	for (int i = 1; i <= m; i++)
	{
		if (eps[i] >= s)
			s = eps[i];
	}
	return s;
}

void grafic(double U[], double t[], int k)
{
	setlocale(LC_NUMERIC, "C");

	FILE *name;
	name = fopen("D:\\new.txt", "w");

	for (int j = 0; j <= N; j++)
	{
		fprintf(name, "%f %f %f\n", x[j], U[j], u(x[j], t[k]));
	}

	fclose(name);

	FILE *name1;
	name1 = fopen("D:\\new1.txt", "w");
	double xs;

	for (int i = 0; i <= K; i++)
	{
		xs = i*l / K;

		fprintf(name1, "%f %f\n", xs, u(xs, t[k]));
	}

	fclose(name1);

}


int main()
{
	setlocale(LC_NUMERIC, "C");
	setlocale(LC_ALL, "rus");

	double M;
	double r = T * 1. / N;
	M = time(r);
	double *t = new double[M + 1];
	double *eps = new double[M + 1];
	nodes(x, N);

	for (int i = 0; i <= M; i++)
	{
		t[i] = i * r;
	}




	for (int j = 1; j <= N - 1; j++)  //задали значение функции на 0 и 1 слое дл€ шага h
	{
		_U[j] = u0(x[j]);
		U[j] = u0(x[j]) + r * (nu*(u0(x[j - 1]) - 2. * u0(x[j]) + u0(x[j + 1])) / pow(l / N, 2.) + f(x[j], t[0]));
	}
	_U[0] = u0(x[0]);
	_U[N] = u0(x[N]);
	U[0] = yO(t[1]);
	U[N] = yl(t[1]);



	//Ќаходим численное решение и погрешность решени€ дл€ шага h
	eps[0] = max(_U, x, t, 0, N);
	eps[1] = max(U, x, t, 1, N);
	double h = x[1] - x[0];
	for (int k = 1; k <= M - 1; k++)
	{
		d[0] = yO(t[k + 1]);
		d[N] = yl(t[k + 1]);
		for (int i = 1; i <= N - 1; i++)
		{
			d[i] = -f(x[i], t[k]) - nu*((1 - 2 * sigma)*(U[i + 1] - 2 * U[i] + U[i - 1]) / pow(h, 2.) + sigma*(_U[i + 1] - 2 * _U[i] + _U[i - 1]) / pow(h, 2.)) - (_U[i]) / (2 * r);
		}
		solution(U_, x, t, d, k, N);
		eps[k + 1] = max(U_, x, t, k + 1, N);
		for (int i = 0; i <= N; i++)
		{
			_U[i] = U[i];
			U[i] = U_[i];
		}
	}



	grafic(U, t, M);
	cout << maxEps(eps, M) << endl;
	_getch();
	return 0;
}

