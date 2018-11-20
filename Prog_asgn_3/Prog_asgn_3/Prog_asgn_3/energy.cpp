//=================================================
// Source file
//=================================================
#include "constant.h"
#include "print_fcns.h"
#include "error_fcns.h"
#include "exact.h"
#include "thomas.h"

void ghost(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int i = 1; i < imax-1; i++)
	{
		T[0][i] = -T[1][i];
		T[jmax-1][i] = 2 - T[jmax-2][i];
	}
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		T[j][0] = 2*(y + 0.75*Pr*Ec*pow(ubar, 2)*(1 - pow((1 - 2*y), 4))) - T[j][1];
		// T[j][imax-1] = 2*(y + 0.75*Pr*Ec*pow(ubar, 2)*(1 - pow((1 - 2*y), 4))) - T[j][imax-2];
		T[j][imax-1] = T[j][imax-2];
	}
}

void init(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			T[j][i] = y;
			u[j][i] = 6.0*ubar*y*(1-y);
			v[j][i] = 0;
		}
	}
	ghost(T, u, v);
}

vector<vector<double>> source(vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> S(jmax, vector<double>(imax));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			S[j][i] = Ec / Re * (2 * pow((u[j][i + 1] - u[j][i - 1]) / (2 * dx), 2) 
								+ 2 * pow((v[j + 1][i] - v[j - 1][i]) / (2 * dy), 2) 
								+ pow((v[j][i + 1] - v[j][i - 1]) / (2 * dx) + (u[j + 1][i] - u[j - 1][i]) / (2 * dy), 2));
		}
	}
	return S;
}

vector<vector<double>> FI2C(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &S)
{
	vector<vector<double>> FI(jmax, vector<double>(imax));
	vector<vector<double>> con(jmax, vector<double>(imax));
	vector<vector<double>> dif(jmax, vector<double>(imax));

	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			// Convective term
			con[j][i] = -(u[j][i + 1] * T[j][i + 1] - u[j][i - 1] * T[j][i - 1]) / (2 * dx) - (v[j + 1][i] * T[j + 1][i] - v[j - 1][i] * T[j - 1][i]) / (2 * dy);

			// Diffusive term
			dif[j][i] = ((T[j][i + 1] - 2 * T[j][i] + T[j][i - 1]) / pow(dx, 2) + (T[j + 1][i] - 2 * T[j][i] + T[j - 1][i]) / pow(dy, 2)) / (Re * Pr);

			FI[j][i] = (con[j][i] + dif[j][i]) + S[j][i];
		}
	}
	return FI;
}

void EE(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	init(T, u, v);
	vector<vector<double>> S = source(u, v);
	vector<vector<double>> T0(jmax, vector<double>(imax));
	int it = 0;
	double delta = 1.0;
	while (abs(delta) > pow(10, -5))
	{
		it++;
		vector<vector<double>> FI = FI2C(T, u, v, S);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T0[j][i] = T[j][i];
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		ghost(T, u, v);

		delta = maxChange(T0, T);
	}
	// printVec2D(T);
	// 
	vector<vector<double>> ExT = exactTemp();
	vector<vector<double>> err = error(T, ExT);
	double L2 = L2Norm(err);

	cout << "Domain: " << imax-2 << " X " << jmax-2 << endl;
	cout << setprecision(6) << "L2 norm: " << L2 << endl;
	cout << "Iterations: " << it << endl; 
}

void RK2(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	init(T, u, v);
	vector<vector<double>> S = source(u, v);
	vector<vector<double>> T0(jmax, vector<double>(imax));
	int it = 0;
	double delta = 1.0;
	while (abs(delta) > pow(10, -5))
	{
		it++;
		// Intermediate Step
		vector<vector<double>> FIint = FI2C(T, u, v, S);
		vector<vector<double>> Tint(jmax, vector<double>(imax));
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T0[j][i] = T[j][i];
				Tint[j][i] = T[j][i] + dt/2.0 * FIint[j][i];
			}
		}
		ghost(Tint, u, v);

		// Full Step
		vector<vector<double>> FI = FI2C(Tint, u, v, S);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] = T[j][i] + dt * FI[j][i];
			}
		}
		ghost(T, u, v);

		delta = maxChange(T0, T);
	}
	// printVec2D(T);
	cout << setprecision(6) << delta << endl;
	vector<vector<double>> ExT = exactTemp();
	vector<vector<double>> err = error(T, ExT);
	double L2 = L2Norm(err);
	cout << setprecision(6) << "L2 norm: " << L2 << endl;
	cout << "Iterations: " << it << endl; 
}

void Imp(vector<vector<double>> &T, vector<vector<double>> &u, vector<vector<double>> &v)
{
	double dt = 0.1;
	init(T, u, v);
	printVec2D(T);
	vector<vector<double>> S = source(u, v);		// Constant
	vector<vector<double>> T0(jmax, vector<double>(imax));
	int it = 0;
	double delta = 1.0;
	while (abs(delta) > pow(10, -5))
	{
		it++;
		vector<vector<double>> FI = FI2C(T, u, v, S);	// Only changes on every time loop

		// Matrix {[I] + dt*[Dx]}
		vector<vector<double>> DX(imax, vector<double>(3));
		vector<double> FIx(imax);
		vector<vector<double>> Ttilda(jmax);
		double ax = dt / (Re*Pr*pow(dx, 2));
		double bx = dt / (2 * dx);

		for (int j = 1; j < jmax - 1; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				T0[j][i] = T[j][i];

				DX[i][0] = -bx * u[j][i-1] - ax;

				DX[i][1] = 1 + 2 * ax;

				DX[i][2] = bx * u[j][i+1] - ax;

				FIx[i] = dt * FI[j][i];
			}

			DX[0][0] = 1;
			DX[0][1] = 1;
			DX[imax-1][1] = 1;
			DX[imax-1][2] = -1;

			FIx[0] = 0;
			FIx[imax-1] = 0;

			SolveThomas(DX, FIx, imax);
			Ttilda[j] = FIx;
		}
		//printVec2D(Ttilda);

		// Matrix {[I] + dt*[Dy]}
		vector<vector<double>> DY(jmax, vector<double>(3));
		vector<double> Tty(jmax);		// Ttilda in vector form for each column
		vector<vector<double>> deltaT(imax);
		double ay = dt / (Re*Pr*pow(dy, 2));
		double by = dt / (2 * dy);

		for (int i = 1; i < imax -1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				DY[j][0] = -by * v[j-1][i] - ay;

				DY[j][1] = 1 + 2 * ay;

				DY[j][2] = by * v[j+1][i] - ay;

				Tty[j] = Ttilda[j][i];
			}

			DY[0][0] = 1;
			DY[0][1] = 1;
			DY[jmax-1][1] = 1;
			DY[jmax-1][2] = -1;

			Tty[0] = 0;
			Tty[jmax-1] = 0;

			SolveThomas(DY, Tty, jmax);
			deltaT[i] = Tty;
		}
		//printVec2D(deltaT);
		
		//printVec2D(T);
		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				T[j][i] += deltaT[i][j];
			}
		}
		ghost(T, u, v);

		delta = maxChange(T0, T);
		
	}
	printVec2D(T);

	vector<vector<double>> ExT = exactTemp();
	cout << setprecision(6) << delta << endl;
	vector<vector<double>> err = error(T, ExT);
	double L2 = L2Norm(err);
	cout << setprecision(6) << "L2 norm: " << L2 << endl;
	cout << "Iterations: " << it << endl; 
}



int main()
{
	vector<vector<double>> T(jmax, vector<double>(imax));
	vector<vector<double>> u(jmax, vector<double>(imax));
	vector<vector<double>> v(jmax, vector<double>(imax));
	


	//init(T, u, v);

	//vector<vector<double>> FI = FI2C(T, u, v, S);
	//printVec2D(FI);

	int start_s=clock();
	// Imp(T, u, v);

	// RK2(T, u, v);

	EE(T, u, v);
	int stop_s=clock();
	cout << "time [sec]: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;


	// the code you wish to time goes here



	// vec2File("T.dat",T);

	//vector<vector<double>> ExT = exactTemp();
	// vector<vector<double>> TE = error(T, ExT);
	// double L2TE = L2Norm(TE);

	//printVec2D(ExT);
	// printVec2D(TE);
	// cout << L2TE;
	
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

