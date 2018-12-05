//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_fcns.h"
#include "error_fcns.h"
#include "exact.h"
#include "thomas.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<double>> &P, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);

	}
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);

	}
}

// Initializes the domain
void init(vector<vector<double>> &P, vector<vector<double>> &u, vector<vector<double>> &v)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			P[j][i] = y;
			u[j][i] = 6.0*ubar*y*(1-y);
			v[j][i] = 0;
		}
	}
	ghost(P, u, v);
}

vector<vector<double>> Pflux(vector<vector<double>> &P, vector<vector<double>> &u, vector<vector<double>> &v)
{
	vector<vector<double>> PFI(jmax, vector<double>(imax));
	vector<vector<double>> F(jmax, vector<double>(imax));
	vector<vector<double>> G(jmax, vector<double>(imax));
	
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			
		}
	}
	return FI;


}


int main()
{
	vector<vector<double>> P(jmax, vector<double>(imax));
	vector<vector<double>> u(jmax, vector<double>(imax));
	vector<vector<double>> v(jmax, vector<double>(imax));



	getchar();
	return 0;
}