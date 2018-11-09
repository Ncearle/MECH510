#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;

constexpr int imax = 81;
constexpr double xmax = 1.0;
constexpr double dx = xmax / (imax - 1);

constexpr double u = 2.0;
constexpr double CFL = 0.5 * 0.8;
constexpr double dt = CFL * dx / u;
//constexpr double dt = 0.0101;		// CFL and dt for part 2
//constexpr double CFL = u * dt / dx;
constexpr double tend = 1.0;
constexpr int tmax = tend / dt + 1;	// Adding 1 to counter the loss when converting double to int

void printVector(vector<double> vec) // Prints a vector to the terminal
{
	for (int i = 0; i < vec.size(); i++)
	{
		cout << vec[i] << "  ";
	}
	cout << "\n\n";
}

void setInit(vector<double> &vec)		// Set initial conditions
{
	for (int i = 0; i < imax; i++)
	{
		vec[i] = -sin(2 * M_PI*(i - 0.5) / (imax - 1));		// Part 1
	}

	//for (int i = 1; i < (imax - 1) / xmax + 1; i++)
	//{
	//	vec[i] = (i - 0.5)*xmax / (imax - 1);				// Part 2
	//}
	
}

vector<double> exact(double t)			// Calculates the exact solution at any time
{
	vector<double> exact(imax, 0.0);
	for (int i = 0; i < imax; i++)
	{
		double x = (i - 0.5) / (imax - 1);
		exact[i] = sin(2 * M_PI*(2 * t - x));
	}
	return exact;
}

double L1Norm(vector<double> err)		// Calculates the L1 Norm of the solution given an error vector
{
	double sum = 0.0;
	for (int i = 1; i < imax; i++)
	{
		sum += abs(err[i]);
	}
	double L1 = sum / (imax - 1);
	return L1;
}

double L2Norm(vector<double> err)		// Calculates the L2 Norm of the solution given an error vector
{
	double sum = 0.0;
	for (int i = 1; i < imax; i++)
	{
		sum += pow(err[i], 2);
	}
	double L2 = sqrt(sum / (imax - 1));
	return L2;
}

double Linf(vector<double> err)		// Calculates the Linfinite norm of the solution given an error vector
{
	double Linf = 0.0;
	for (int i = 1; i < imax; i++)
	{
		if (abs(err[i]) > abs(Linf))
		{
			Linf = abs(err[i]);
		}
	}
	return Linf;
}

vector<double> FI1U(vector<double> &vec, double t)	// 1st order upwind flux integral
{
	vector<double> temp(imax, 0.0);
	for (int i = 1; i < vec.size(); i++)
	{
		temp[i] = -u*(vec[i] - vec[i - 1]) / dx;
	}
	return temp;
}

vector<double> FI2U(vector<double> &vec, double t) // 2nd order upwind flux integral
{
	vector<double> temp(imax, 0.0);
	for (int i = 2; i < vec.size(); i++)
	{
		temp[i] = -u*(3 * vec[i] - 4 * vec[i - 1] + vec[i - 2]) / (2 * dx);
	}
	return temp;
}

void setBound1(vector<double> &vec, double t)	// Boundary Conditions for first order scheme
{

	vec[0] = 2 * sin(4 * M_PI * t) - vec[1];
}

void setBound2(vector<double> &vec1, vector<double> &vec0, double t)		// Boundary conditions for second order scheme
{
		double t32 = (3.0 * vec0[1] - vec0[0]) / 2.0;						// Flux across i_3/2 for the step
		double t12 = sin(4.0 * M_PI*(t - dt / 2.0));						// Flux across i_1/2 for the step

		vec1[1] = vec0[1] - u * dt / 2.0 * (t32 - t12) / dx;				// Set first cell value using fluxes above
		vec1[0] = 2.0 * sin(4.0 * M_PI * (t)) - vec1[1];					// Set ghost cell using new first cell value
}

void EE(vector<double> &vec, vector<double> &err, vector<double> &L2, vector<double> &L1, vector<double> &LI, int order) // Explicit Euler time advance
{
	setInit(vec);

	// Stuff to write data to csv files
	/*ofstream data;
	ofstream error;
	ofstream L2N;
	ofstream L1N;
	ofstream LIN;
	data.open("Data_EE_FI1U_80.csv");
	error.open("Error_EE_FI1U_80.csv");
	L2N.open("L2Norm_EE_FI1u_80.csv");
	L1N.open("L1Norm_EE_FI1u_80.csv");
	LIN.open("LInf_EE_FI1U_80.csv");*/
	
	for (int n = 1; n < tmax + 1; n++)
	{
		double t = n * dt;
		if (order == 1)				// First order upwind scheme for EE
		{	
			vector<double> temp = FI1U(vec, t);
			for (int i = 1; i < imax; i++)
			{
				vec[i] = vec[i] + dt * temp[i];
			}
			setBound1(vec, t);
		}
		else if (order == 2)			// Second order upwind scheme for EE
		{
			vector<double> temp = FI2U(vec, t);
			for (int i = 2; i < imax; i++)
			{
				vec[i] = vec[i] + dt * temp[i];
			}
			setBound2(vec, vec, t);
		}
		else
		{
			cout << "ERROR!!! Order must be 1 or 2";
			break;
		}

		// Exact Solution and Error Calculation
		vector<double> T_exact = exact(t);
		for (int i = 1; i < imax; i++)
		{
			err[i] = T_exact[i] - vec[i];	// Calculates error as difference between the exact solution and numerical solution

			//data << vec[i] << ", ";		// Print data to a csv file
			//error << err[i] << ", ";		// Print error to a csv file
		}
		L2[n - 1] = L2Norm(err);
		L1[n - 1] = L1Norm(err);
		LI[n - 1] = Linf(err);

		/*data << endl;
		error << endl;
		L2N << L2[n - 1] << endl;
		L1N << L1[n - 1] << endl;
		LIN << LI[n - 1] << endl;*/

	}
	printVector(vec);
	/*data.close();
	error.close();
	L2N.close();
	L1N.close();
	LIN.close();*/
}

void RK2(vector<double> &vec, vector<double> &err, vector<double> &L2, vector<double> &L1, vector<double> &LI, int order)		// Runge Kutta 2 time advance
{
	setInit(vec);

	// Stuff to write data to csv files
	/*ofstream data;
	ofstream error;
	ofstream L2N;
	ofstream L1N;
	ofstream LIN;
	data.open("Data_RK2_FI2U_80.csv");
	error.open("Error_RK2_FI2U_80.csv");
	L2N.open("L2Norm_RK2_FI2U_80.csv");
	L1N.open("L1Norm_RK2_FI2U_80.csv");
	LIN.open("LInf_RK2_FI2U_80.csv");*/

	for (int n = 1; n < tmax+1; n++)
	{
		double t = n * dt;
		if (order == 1)				// Firsts order Upwind Scheme with RK2
		{
			// Intermediate step
			vector<double> temp1 = FI1U(vec, t - dt / 2.0);
			vector<double> w_int(imax, 0.0);
			for (int i = 1; i < imax; i++)
			{
				w_int[i] = vec[i] + (dt / 2.0) * temp1[i];
			}
			setBound1(w_int, t - dt / 2);

			// Full step
			vector<double> temp2 = FI1U(w_int, t);
			for (int i = 1; i < imax; i++)
			{
				vec[i] = vec[i] + dt * temp2[i];
			}
			setBound1(vec, t);
		}
		else if (order == 2)			// Second order Upwind Scheme with RK2
		{
			// Intermediate step
			vector<double> temp1 = FI2U(vec, t - dt / 2.0);
			vector<double> w_int(imax, 0.0);
			for (int i = 2; i < imax; i++)
			{
				w_int[i] = vec[i] + (dt / 2.0) * temp1[i];
			}
			setBound2(w_int, vec, t - dt / 2);

			// Full step
			vector<double> temp2 = FI2U(w_int, t);
			for (int i = 2; i < imax; i++)
			{
				vec[i] = vec[i] + dt * temp2[i];
			}
			setBound2(vec, w_int, t);
		}
		else
		{
			cout << "ERROR!!! Order must be 1 or 2";
			break;
		}
		// Exact Solution and error
		vector<double> T_exact = exact(t);

		for (int i = 1; i < imax; i++)
		{
			err[i] = T_exact[i] - vec[i];	// Calculates error as difference between the exact solution and numerical solution

			//data << vec[i] << ", ";		// Print data to a csv file
			//error << err[i] << ", ";		// Print error to a csv file
		}
		L2[n - 1] = L2Norm(err);
		L1[n - 1] = L1Norm(err);
		LI[n - 1] = Linf(err);
		
		/*data << endl;
		error << endl;
		L2N << L2[n - 1] << endl;
		L1N << L1[n - 1] << endl;
		LIN << LI[n - 1] << endl;*/
	   
	}
	/*data.close();
	error.close();
	L2N.close();
	L1N.close();
	LIN.close();*/
}

int main()
{

	vector<double> T(imax,0.0);			// Temperature vector, for each cell
	vector<double> E(imax, 0.0);			// Error vector, for each cell
	vector<double> L2(tmax, 0.0);			// L2 Norm vector, for each timestep
	vector<double> L1(tmax, 0.0);			// L1 Norm vector
	vector<double> LI(tmax, 0.0);			// L infinite norm
	vector<double> T_exact = exact(tend);	// Exact soltution at given timestep
	
	RK2(T, E, L2, L1, LI, 2);

	//EE(T, E, L2, L1, LI, 1);

	//printVector(T);
	//printVector(T_exact);
	
	cout << "\n\nL1Norm: " << L1[tmax-1] << "\n\n";
	cout << "\n\nL2Norm: " << L2[tmax - 1] << "\n\n";
	cout << "\n\nLInfinite: " << LI[tmax - 1] << "\n\n";

	//cout << "\n\nTime-step: " << dt << endl;
	//cout << "\nCFL#: " << CFL << endl;

	getchar();
	return 0;
}