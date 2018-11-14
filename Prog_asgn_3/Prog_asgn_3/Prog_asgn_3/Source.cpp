#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include "energy.h"

using namespace std;

// Domain Constants
constexpr int imax = 12;
constexpr int jmax = 12;
constexpr double xmax = 1;
constexpr double ymax = 1;
constexpr double dx = xmax / (imax - 2);
constexpr double dy = ymax / (jmax - 2);


// Constants
constexpr double pi = M_PI;	// Pi
constexpr int Re = 50;		// Reynolds number
constexpr double Pr = 0.7;	// Prandtl number
constexpr double Ec = 0.1;	// Eckert number

constexpr int ubar = 3;		// average velocity in x [m/s]

constexpr int u0 = 1;		// Initial x velocity
constexpr int v0 = 1;		// Initial y velocity
constexpr int T0 = 1;		// Initial temperature

void init(double(&T)[jmax][imax], double(&u)[jmax][imax], double(&v)[jmax][imax])
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			
			T[j][i] = T0 * cos(pi*x) * sin(pi*y);
			u[j][i] = u0 * y * sin(pi*x);
			v[j][i] = v0 * x * cos(pi*y);
		}
	}
}

void FI2C(double(&T)[jmax][imax], double(&u)[jmax][imax], double(&v)[jmax][imax])
{
	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			// Convective term
			double C = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - u[j - 1][i] * T[j - 1][i]) / (2 * dx);

			// Diffusive term
			double D = 1 / (Re*Pr) * ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2));

			temp[j][i] = C + D;
		}
	}
}

void exact(double(&T)[jmax][imax])
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			T[j][i] = u0 * T0*pi*cos(2 * pi*x)*y*sin(pi*y) + v0 * T0*pi*x*cos(pi*x)*cos(2 * pi*y) + 2 * T0*pow(pi, 2)*cos(pi*x)*sin(pi*y) / (Re*Pr);
		}
	}
}


int main()
{
	double u[jmax][imax] = { {0} };
	double v[jmax][imax] = { {0} };
	double T[jmax][imax] = { {0} };
	double Exact[jmax][imax] = { {0} };

	init(T, u, v);
	printMatrix(T, jmax, imax, "T");
	printMatrix(u, jmax, imax, "u");
	printMatrix(v, jmax, imax, "v");
	double FI = FI2C(T, u, v);
	exact(Exact);

	printMatrix(FI, jmax, imax, "FI");
	printMatrix(Exact, jmax, imax, "Exact");
	

	getchar();
	return 0;	
}

