#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include "print_fcns.h"
#include "error_fcns.h"

using namespace std;

// Domain Constants

constexpr double xmax = 1;
constexpr double ymax = 1;
constexpr double dx = xmax / (imax - 2);
constexpr double dy = ymax / (jmax - 2);

constexpr int ubar = 3;		// average velocity in x [m/s]
constexpr double CFL = 0.4;	// CFL 80% of maximum for RK2
constexpr double dt = CFL * dx / ubar;	//Time step

constexpr double tend = 1.0;
constexpr int tmax = tend / dt + 1; //Adding 1 to counter the loss when converting to int

// Constants
constexpr double pi = M_PI;	// Pi
constexpr int Re = 50;		// Reynolds number
constexpr double Pr = 0.7;	// Prandtl number
constexpr double Ec = 0.1;	// Eckert number

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

double** source(double(&u)[jmax][imax], double(&v)[jmax][imax])
{
	double** S = 0;
	S = new double*[jmax];
	for (int j = 0; j < jmax; j++)
	{
		S[j] = new double[imax];
		for (int i = 0; i < imax; i++)
		{

			S[j][i] = Ec/Re * (2*pow((u[j][i+1] - u[j][i-1])/(2*dx),2) + 2*pow((v[j+1][i] - v[j-1][i])/(2*dy),2) + pow((v[j][i+1] - v[j][i-1])/(2*dx) + (u[j+1][i] - u[j-1][i])/(2*dy),2));
		}
	}
	return S;
}

double** FI2C(double(&T)[jmax][imax], double(&u)[jmax][imax], double(&v)[jmax][imax])
{
	double** temp = 0;
	temp = new double*[jmax];
	double Con[jmax][imax] = {{0}};
	double Dif[jmax][imax] = {{0}};
	double** S = source(u, v);
	for (int j = 1; j < jmax - 1; j++)
	{
		temp[j] = new double[imax];
		for (int i = 1; i < imax - 1; i++)
		{
			// Convective term
			Con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dx);

			// Diffusive term
			Dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			temp[j][i] = -(Con[j][i] + Dif[j][i]) + S[i][j];
		}
	}
	return temp;
}

double** exactFlux()
{
	double** EF = 0;
	EF = new double*[jmax];
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		EF[j] = new double[imax];
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			EF[j][i] = u0 * T0*pi*cos(2 * pi*x)*y*sin(pi*y) + v0 * T0*pi*x*cos(pi*x)*cos(2 * pi*y) + 2 * T0*pow(pi, 2)*cos(pi*x)*sin(pi*y) / (Re*Pr);
		}
	}
	return EF;
}

double** exactSource()
{
	double** ES = 0;
	ES = new double*[jmax];
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		ES[j] = new double[imax];
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			ES[j][i] = Ec/Re * (2*pow(u0*pi*cos(pi*x)*y,2) + 2*pow(v0*pi*x*sin(pi*y),2) + pow(u0*sin(pi*x) + v0*cos(pi*y),2));
		}
	}
	return ES;
}

void RK2(double(&T)[jmax][imax], double(&u)[jmax][imax], double(&v)[jmax][imax])
{
	init(T, u, v);

	for (int n = 1; n < tmax + 1; n++)
	{
		double t = n * dt;
		// Intermediate Step
		double** FIint = FI2C(T, u, v);
		double Tint[jmax][imax] = { {0} };
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				Tint[j][i] = T[j][i] + dt/2 * FIint[j][i];
			}
		}

		// Full Step
		double** FI = FI2C(Tint, u, v);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
	}
}

int main()
{
	double u[jmax][imax] = { {0} };
	double v[jmax][imax] = { {0} };
	double T[jmax][imax] = { {0} };

	RK2(T, u, v);

	printMatrix(T, jmax, imax, "T");

	// init(T, u, v);

	// double** FI =  FI2C(T, u, v);
	// double** S = source(u, v);

	// double** EFI = exactFlux();
	// double** ES = exactSource();

	// double** fluxError = error(FI, EFI);
	// double** sourceError = error(S, ES);

	// double fluxL2 = L2Norm(fluxError);
	// double sourceL2 = L2Norm(sourceError);
	
	// printMatrix(FI, jmax, imax, "Flux Integral");

	// printMatrix(EFI, jmax, imax, "Exact Flux Term");
	// // printMatrix(S, jmax, imax, "Source Term");
	// // printMatrix(ES, jmax, imax, "Exact Source Term");

	// printMatrix(fluxError, jmax, imax, "Flux Integral error");
	// cout << "\nFlux L2Norm = " << fluxL2 << endl;

	// // printMatrix(sourceError, jmax, imax, "Source Term error");
	// // cout << "\nSource L2Norm = " << sourceL2 << endl;

	getchar();
	return 0;
}

