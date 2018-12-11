//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "thomas.h"
#include "vecmat_funcs.h"

// Given domain and velocities, sets and updates the bounary conditions
void ghost(vector<vector<vector<double>>> &U)
{
	// Top and Bottom Ghost Cells
	for (int i = 1; i < imax-1; i++)
	{
		double x = (i - 0.5) / (imax - 2);
		// pressure -- Neumann boundary condition (no flux across boundary)
		U[0][i][0] = U[1][i][0];
		U[jmax-1][i][0] = U[jmax-2][i][0];

		// u velocity -- Dirichlet boundary condition (no-slip at wall)
		U[0][i][1] = -U[1][i][1];
		U[jmax-1][i][1] = -U[jmax-2][i][1];

		// v velocity -- Dirichlet boundary condition (no-slip at wall)
		U[0][i][2] = -U[1][i][2];
		U[jmax-1][i][2] = -U[jmax-2][i][2];
	}

	// Left and Right Ghost Cells
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		// pressure -- Neumann boundary condition (no flux across boundary)
		U[j][0][0] = U[j][1][0];		
		U[j][imax-1][0] = U[j][imax-2][0];

		// u velocity -- Dirichlet boundary condition (no-slip at wall)
		U[j][0][1] = -U[j][1][1];
		U[j][imax-1][1] = -U[j][imax-2][1];

		// v velocity -- Dirichlet boundary condition (no-slip at wall)
		U[j][0][2] = -U[j][1][2];
		U[j][imax-1][2] = -U[j][imax-2][2];
	}
}

// Initializes the domain (Pressure and Velocity)
void init(vector<vector<vector<double>>> &U)
{
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);

			// Pressure
			U[j][i][0] = P0 * cos(pi*x) * cos(pi * y);

			// u velocity
			U[j][i][1] = u0 * sin(pi*x) * sin(2*pi*y);

			// v velocity
			U[j][i][2] = v0 * sin(2*pi*x) * sin(pi*y);
		}
	}
	ghost(U);
}

vector<vector<vector<double>>> flux(vector<vector<vector<double>>> &U)
{
	vector<vector<vector<double>>> flux(jmax, vector<vector<double>>(imax, vector<double>(3)));
	for (int j = 1; j < jmax-1; j++)
	{
		for (int i = 1; i < imax-1; i++)
		{
			vector<double> Fp(3);	// F_i+1/2,j
			vector<double> Fm(3);	// F_i-1/2,j
			vector<double> Gp(3);	// G_i,j+1/2
			vector<double> Gm(3);	// G_i,j-1/2

			Fp[0] = (U[j][i+1][1] + U[j][i][1]) / (2*B);
			Fp[1] = pow((U[j][i+1][1] + U[j][i][1])/2, 2) + (U[j][i+1][0] + U[j][i][0])/2 - (U[j][i+1][1] - U[j][i][1])/(Re*dx);
			Fp[2] = ((U[j][i+1][1] + U[j][i][1]) / 2) * ((U[j][i+1][2] + U[j][i][2]) / 2) - (U[j][i+1][2] - U[j][i][2])/(Re*dx);

			Fm[0] = (U[j][i][1] + U[j][i-1][1]) / (2*B);
			Fm[1] = pow((U[j][i][1] + U[j][i-1][1])/2, 2) + (U[j][i][0] + U[j][i-1][0])/2 - (U[j][i][1] - U[j][i-1][1])/(Re*dx);
			Fm[2] = ((U[j][i][1] + U[j][i-1][1]) / 2) * ((U[j][i][2] + U[j][i-1][2]) / 2) - (U[j][i][2] - U[j][i-1][2])/(Re*dx);

			Gp[0] = (U[j+1][i][2] + U[j][i][2]) / (2*B);
			Gp[1] = ((U[j+1][i][1] + U[j][i][1]) / 2) * ((U[j+1][i][2] + U[j][i][2]) / 2) - (U[j+1][i][1] - U[j][i][1])/(Re*dy);
			Gp[2] = pow((U[j+1][i][2] + U[j][i][2])/2, 2) + (U[j+1][i][0] + U[j][i][0])/2 - (U[j+1][i][2] - U[j][i][2])/(Re*dy);

			Gm[0] = (U[j][i][2] + U[j-1][i][2]) / (2*B);
			Gm[1] = ((U[j][i][1] + U[j-1][i][1]) / 2) * ((U[j][i][2] + U[j-1][i][2]) / 2) - (U[j][i][1] - U[j-1][i][1])/(Re*dy);
			Gm[2] = pow((U[j][i][2] + U[j-1][i][2])/2, 2) + (U[j][i][0] + U[j-1][i][0])/2 - (U[j][i][2] - U[j-1][i][2])/(Re*dy);


			for (int k = 0; k < 3; k++)
			{
				flux[j][i][k] = - (Fp[k] - Fm[k])/dx - (Gp[k] - Gm[k])/dy;
			}
		}
	}
	return flux;
}

vector<vector<double>> jac(vector<vector<vector<double>>> &U, string FG, int FGpm, int Upm, int j, int i)
{
	vector<vector<double>> J(3, vector<double>(3));
	if (FG == "F")
	{
		J[0][1] = -1.0 / (2*B);
		J[1][0] = -1.0 / 2.0;
		J[1][1] = -(U[j][i][1] + U[j][i+FGpm][1])/2 + Upm / (Re*dx);
		J[2][1] = -(U[j][i][2] + U[j][i+FGpm][2])/4;
		J[2][2] = -(U[j][i][1] + U[j][i+FGpm][1])/4 + Upm / (Re*dx);
	}

	else if (FG == "G")
	{
		J[0][2] = -1.0 / (2*B);
		J[1][1] = -(U[j][i][2] + U[j+FGpm][i][2])/4 + Upm / (Re*dy);
		J[1][2] = -(U[j][i][1] + U[j+FGpm][i][1])/4;
		J[2][0] = -1.0 / 2.0;
		J[2][2] = -(U[j][i][2] + U[j+FGpm][i][2])/2 + Upm / (Re*dy);
	}
	return J;
}

void Imp(vector<vector<vector<double>>> &U)
{


}

int main()
{
	vector<vector<vector<double>>> U(jmax, vector<vector<double>>(imax, vector<double>(3)));
	init(U);

	vector<vector<vector<double>>> FI = flux(U);
	vector<vector<vector<double>>> ExFlux = exactFlux();
	vector<vector<vector<vector<double>>>> DX(imax, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(3))));
	vector<vector<double>> FIx(imax, vector<double>(3));

	int j = 10;
	int i = 10;

	vector<vector<vector<double>>> Udelta(jmax, vector<vector<double>>(imax, vector<double>(3)));
	Udelta[j][i][0] = pow(10, -6);
	Udelta[j][i][1] = pow(10, -6);
	Udelta[j][i][2] = pow(10, -6);

	vector<vector<vector<double>>> Unew = Madd3D(U, Udelta);
	vector<vector<vector<double>>> FInew = flux(Unew);

	vector<vector<double>> Ax = jac(U, "F", -1, -1, j, i);
	vector<vector<double>> Bxp = jac(U, "F", 1, -1, j, i);
	vector<vector<double>> Bxm = jac(U, "F", -1, 1, j, i);
	vector<vector<double>> Cx = jac(U, "F", 1, 1, j, i);

	vector<vector<double>> Ay = jac(U, "G", -1, -1, j, i);
	vector<vector<double>> Byp = jac(U, "G", 1, -1, j, i);
	vector<vector<double>> Bym = jac(U, "G", -1, 1, j, i);
	vector<vector<double>> Cy = jac(U, "G", 1, 1, j, i);

	Ax = ScaM(-1, Ax);
	Ay = ScaM(-1, Ay);
	vector<vector<double>> Bx = Msub(Bxp, Bxm);
	vector<vector<double>> By = Msub(Byp, Bym);
	Cx = ScaM(1, Cx);
	Cy = ScaM(1, Cy);

	vector<vector<vector<double>>> FIdelta = error(FI, FInew);

	// printVec(FIdelta[j][i]);
	// printVec3D(FI, 0);
	vector<double> X = Vsub(FInew[j][i], FI[j][i]);
	vector<double> Y = MVM(Bxp, Udelta[10][10]);

	Y = ScaV(-1/dx, Y);	

	printVec(X);
	// printVec(Y);

	vector<double> AxU = MVM(Ax, Udelta[j][i-1]);
	vector<double> AyU = MVM(Ay, Udelta[j-1][i]);
	vector<double> BxU = MVM(Bx, Udelta[j][i]);
	vector<double> ByU = MVM(By, Udelta[j][i]);
	vector<double> CxU = MVM(Cx, Udelta[j][i+1]);
	vector<double> CyU = MVM(Cy, Udelta[j+1][i]);
	vector<double> RHS(3);
	for (int n = 0; n < 3; n++)
	{
		RHS[n] = (AxU[n]  + BxU[n] + + CxU[n])/dx + (AyU[n] + ByU[n] + CyU[n])/dy;
	}
	printVec(RHS);

	// printVec(AxU);
	// printVec(AyU);
	// printVec(BxU);
	// printVec(ByU);
	// printVec(CxU);
	// printVec(CyU);

	printVec(Vsub(X, RHS));

	int k = 1;	// Index for solution: 0 = pressure; 1 = u velocity; 2 = v velocity
	// printVec3D(FIdelta, k, p);

	
	// printVec3D(U, k, p);
	// printVec3D(ExFlux, k, p);
	// printVec3D(FI, k, p);

	// vector<vector<vector<double>>> E = error(FI, ExFlux);
	// vector<double> L2 = L2Norm(E);

	// printVec(L2);

	// }

			
	getchar();
	return 0;
}
