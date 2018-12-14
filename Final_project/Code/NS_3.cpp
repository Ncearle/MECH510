//=================================================
// Source file for Navier-Stokes Equation
//=================================================
#include "constant.h"
#include "print_funcs.h"
#include "error_funcs.h"
#include "exact.h"
#include "vecmat_funcs.h"
#include "blocktri_vector.h"

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
	for (int j = 1; j < jmax-1; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 1; i < imax-1; i++)
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
		J[0][1] = 1.0 / (2*B);
		J[1][0] = 1.0 / 2.0;
		J[1][1] = (U[j][i][1] + U[j][i+FGpm][1])/2 + Upm / (Re*dx);
		J[2][1] = (U[j][i][2] + U[j][i+FGpm][2])/4;
		J[2][2] = (U[j][i][1] + U[j][i+FGpm][1])/4 + Upm / (Re*dx);

		J = ScaM(1/dx, J);
	}

	else if (FG == "G")
	{
		J[0][2] = 1.0 / (2*B);
		J[1][1] = (U[j][i][2] + U[j+FGpm][i][2])/4 + Upm / (Re*dy);
		J[1][2] = (U[j][i][1] + U[j+FGpm][i][1])/4;
		J[2][0] = 1.0 / 2.0;
		J[2][2] = (U[j][i][2] + U[j+FGpm][i][2])/2 + Upm / (Re*dy);

		J = ScaM(1/dy, J);
	}
	return J;
}

void Imp(vector<vector<vector<double>>> &U)
{
	init(U);
	vector<vector<double>> I = Id(3);

	for (int n = 0; n < 200; n++)
	{
		vector<vector<vector<double>>> FI = flux(U);
		// cout << "Flux Integral: " << endl;
		// printTable(FI);

		vector<vector<vector<vector<double>>>> DX(imax, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(3))));
		vector<vector<double>> FIx(imax, vector<double>(3));
		vector<vector<vector<double>>> Utilda(jmax, vector<vector<double>>(imax, vector<double>(3)));

		for (int j = 1; j < jmax - 1; j++)
		{
			for (int i = 1; i < imax - 1; i++)
			{
				vector<vector<double>> Ax = jac(U, "F", -1, -1, j, i);
				vector<vector<double>> Bxp = jac(U, "F", 1, -1, j, i);
				vector<vector<double>> Bxm = jac(U, "F", -1, 1, j, i);
				vector<vector<double>> Cx = jac(U, "F", 1, 1, j, i);
				vector<vector<double>> Bx = Msub(Bxp, Bxm);

				Ax = ScaM(-dt, Ax);
				Bx = ScaM(dt, Bx);
				Cx = ScaM(dt, Cx);
				
				DX[i][0] = Ax;
				DX[i][1] = Madd(I, Bx);
				DX[i][2] = Cx;

				FIx[i] = ScaV(dt, FI[j][i]);
			}

			DX[0][0] = I;
			DX[0][1] = I;
			DX[0][1][0][0] = -1;

			DX[imax-1][1] = I;
			DX[imax-1][2] = I;
			DX[imax-1][2][0][0] = -1;

			// FIx boundary cells are already set to zero
			SolveBlockTri(DX, FIx, imax);
			Utilda[j] = FIx;
		}
		// cout << "Solving for lines of constant J:" << endl;
		// printTable(Utilda);

		vector<vector<vector<vector<double>>>> DY(jmax, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(3))));
		vector<vector<double>> Uty(jmax, vector<double>(3));
		vector<vector<vector<double>>> deltaU(imax, vector<vector<double>>(jmax, vector<double>(3)));

		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				vector<vector<double>> Ay = jac(U, "G", -1, -1, j, i);
				vector<vector<double>> Byp = jac(U, "G", 1, -1, j, i);
				vector<vector<double>> Bym = jac(U, "G", -1, 1, j, i);
				vector<vector<double>> Cy = jac(U, "G", 1, 1, j, i);
				vector<vector<double>> By = Msub(Byp, Bym);

				Ay = ScaM(-dt, Ay);
				By = ScaM(dt, By);
				Cy = ScaM(dt, Cy);

				DY[j][0] = Ay;
				DY[j][1] = Madd(I, By);
				DY[j][2] = Cy;

				Uty[j] = Utilda[j][i];
			}

			DY[0][1] = I;
			DY[0][2] = I;
			DY[0][2][0][0] = -1;

			DY[jmax-1][0] = I;
			DY[jmax-1][1] = I;
			DY[jmax-1][1][0][0] = -1;

			SolveBlockTri(DY, Uty, jmax);
			deltaU[i] = Uty;
		}
		// cout << "Solving for lines of constant I: " << endl;
		// printTable(deltaU);

		for (int j = 1; j < jmax-1; j++)
		{
			for (int i = 1; i < imax-1; i++)
			{
				for (int k = 0; k < 3; k++)
				{
					U[j][i][k] += deltaU[i][j][k];	// deltaU is in the shape of the transpose of U
				}
			}
		}
		ghost(U);
	}

	printTable(U);
}

int main()
{
	vector<vector<vector<double>>> U(jmax, vector<vector<double>>(imax, vector<double>(3)));
	init(U);
	// cout << "Initial condition: " << endl;
	// printTable(U);
	Imp(U);
	// printTable(U);



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
