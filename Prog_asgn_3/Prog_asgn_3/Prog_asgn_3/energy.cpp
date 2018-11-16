//=================================================
//Source file
//=================================================
#include "constant.h"
#include "print_fcns.h"
#include "error_fcns.h"
#include "exact.h"

void init(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax-1; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			// T[j][i] = T0 * cos(pi*x) * sin(pi*y);	// IC for Part 1,2
			// u[j][i] = u0 * y * sin(pi*x);
			// v[j][i] = v0 * x * cos(pi*y);

			T[j][i] = y;
			u[j][i] = 6*ubar*y*(1-y);
			v[j][i] = 0;
		}
	}
}

void setBound(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int i = 1; i < imax-1; i++)
	{
		T[0][i] = -T[1][i];
		T[jmax-1][i] = 2 - T[jmax-2][i];
	}
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		T[j][imax-1] = 2*(y + 0.75*Pr*Ec*pow(ubar, 2)*(1 - pow((1 - 2*y), 4))) - T[j][imax-2];
		T[j][0] = 2*(y + 0.75*Pr*Ec*pow(ubar, 2)*(1 - pow((1 - 2*y), 4))) - T[j][1];
	}
}

vector<vector<double>> source(vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> S(jmax, vector<double>(imax));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{

			S[j][i] = Ec / Re * (2 * pow((u[j][i + 1] - u[j][i - 1]) / (2 * dx), 2) + 2 * pow((v[j + 1][i] - v[j - 1][i]) / (2 * dy), 2) + pow((v[j][i + 1] - v[j][i - 1]) / (2 * dx) + (u[j + 1][i] - u[j - 1][i]) / (2 * dy), 2));
		}
	}
	return S;
}

vector<vector<double>> FI2C(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> FI(jmax, vector<double>(imax));
	vector<vector<double>> con(jmax, vector<double>(imax));
	vector<vector<double>> dif(jmax, vector<double>(imax));

	vector<vector<double>> S = source(u, v);
	for (int j = 1; j < jmax - 1; j++)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			// Convective term
			con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dx);

			// Diffusive term
			dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			FI[j][i] = (con[j][i] + dif[j][i]) + S[i][j];
		}
	}
	return FI;
}

void EE(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	init(T, u, v);
	setBound(T, u, v);
	printVec2D(T);
	for (int n = 1; n < tmax; n++)
	{
		double t = n * dt;
		vector<vector<double>> FI = FI2C(T, u, v);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		setBound(T, u, v);
		// printVec2D(T);
	}
}

void RK2(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	init(T, u, v);
	setBound(T, u, v);
	printVec2D(T);
	for (int n = 1; n < 2*tmax; n++)
	{
		double t = n * dt;

		// Intermediate Step
		vector<vector<double>> FIint = FI2C(T, u, v);
		vector<vector<double>> Tint(jmax, vector<double>(imax));
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				Tint[j][i] = T[j][i] + dt/2 * FIint[j][i];
			}
		}
		setBound(Tint, u, v);

		// Full Step
		vector<vector<double>> FI = FI2C(Tint, u, v);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		setBound(T, u, v);
		// printVec2D(T);
	}
}

int main()
{
	vector<vector<double>> T(jmax, vector<double>(imax));
	vector<vector<double>> u(jmax, vector<double>(imax));
	vector<vector<double>> v(jmax, vector<double>(imax));

	// init(T, u, v);
	// setBound(T, u, v);


	RK2(T, u, v);
	printVec2D(T);
	// EE(T, u, v);
	// printVec2D(T);

	vec2File("T.dat",T);


	// vector<vector<double>> S = source(u, v);
	// vector<vector<double>> ExS = exactSource();
	// vector<vector<double>> ExFI = exactFlux();
	// vector<vector<double>> FI = FI2C(T, u, v);

	// vector<vector<double>> FE = error(FI, ExFI);
	// vector<vector<double>> SE = error(S, ExS);

	// printVec2D(FE);
	// printVec2D(SE);

	// double L2FE = L2Norm(FE);
	// double L2SE = L2Norm(SE);


	getchar();
	return 0;
}

