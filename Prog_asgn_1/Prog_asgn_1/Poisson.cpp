#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;

constexpr int imax = 102, jmax = 102;
constexpr double xmax = 1.0, ymax = 1.0;
constexpr double dx = xmax / (imax - 2), dy = ymax / (jmax - 2);

// Exact Solution over the domain
double exactSol(double(&exact)[imax][jmax]) {		
	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			exact[j][i] = (cos(M_PI*(i - .5) * xmax/ (imax - 2)) * sinh(M_PI*(j - .5) * ymax/ (jmax - 2))) / sinh(M_PI);
		}
	}
	return **exact;
}

// Calculates the error given the exact and numerical solutions
double error(double (&error)[imax][imax], double(&domain)[imax][jmax], double(&exact)[imax][jmax]) {

	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			error[i][j] = exact[i][j] - domain[i][j];
		}
	}
	return **error;
}

// Calculates the L2 Norm given the error matrix
double L2Norm(double(&error)[imax][jmax]) {			

	double sum = 0;
	
	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			sum += pow(error[i][j], 2);
		}
	}

	return sqrt(sum / ((imax - 1)*(jmax - 2)));
}

// Calculates the maximum change of all cells between iterations given the change matrix
double maxChange(double(&delta)[imax-2][jmax-2]) {			

	double maxChange = 0;

	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			if (abs(delta[i][j]) > maxChange) {
				maxChange = delta[i][j];
			}
		}
	}
	return maxChange;
}

// Boundary conditions for the temperature test case
void setBoundaryT(double(&domain)[imax][jmax]) {
	// Top and Bottom//
	for (int j = 0; j < jmax; j++)
	{
		domain[imax - 1][j] = -domain[imax - 2][j] + 2 * cos(M_PI*(j - 0.5) * ymax / (jmax - 2));
		domain[0][j] = -domain[1][j];
	}
	//Left and Right//
	for (int i = 0; i < imax; i++)
	{
		domain[i][jmax - 1] = domain[i][jmax - 2];
		domain[i][0] = domain[i][1];
	}
}

// Boundary conditions for the Poisson pressure example
void setBoundaryP(double(&domain)[imax][jmax]) {
	// Top and Bottom//
	for (int j = 0; j < jmax; j++)
	{
		domain[imax - 1][j] = -domain[imax - 2][j] + 2 * (5 - 0.5*pow((1 + pow((j - 0.5)*dx, 2)), 3));
		domain[0][j] = domain[1][j];
	}
	//Left and Right//
	for (int i = 0; i < imax; i++)
	{
		domain[i][jmax - 1] = -domain[i][jmax - 2] + 2 * (5 - 0.5*pow((1 + pow((i - 0.5)*dy, 2)), 3));
		domain[i][0] = domain[i][1];
	}
}

// Source term for the Poisson pressure example
void source(double(&S)[imax][jmax]) {
	for (int i = 1; i < imax - 1; i++)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			S[j][i] = -18 * pow(pow((j - 0.5)*dx, 2) + pow((i - 0.5)*dy, 2), 2);
		}
	}
}

// Point Gauss Seidel with over-relaxation
// Runs given a solution matrix (domain), a delta matrix for changes between iterations, the tolerance, a vector to hold the maximum changes at each iteration, the over relaxation factor omega, and the source matrix.
void PGS_OR(double(&domain)[imax][jmax], double(&delta)[imax - 2][jmax - 2], double tol, vector<double> vec, double omega, double(&source)[imax][jmax]){
	double maxDelta = 1.0;
	int it = 1;
	//setBoundaryT(domain);
	setBoundaryP(domain);
	while (maxDelta > tol) {
		for (int i = 1; i < imax - 1; i++)
		{
			for (int j = 1; j < jmax - 1; j++)
			{
				double domain_ij0 = domain[i][j];
				double ddomain_ij = pow(dy,2)/(2*(pow(dx,2)+pow(dy,2))) * (domain[i + 1][j] + domain[i - 1][j]) + pow(dx, 2) / (2 * (pow(dx, 2) + pow(dy, 2)))* (domain[i][j + 1] + domain[i][j - 1]) - (pow(dy,2) * pow(dx, 2) / (2 * (pow(dx, 2) + pow(dy, 2)))) * source[i][j] - domain_ij0;

				domain[i][j] = domain_ij0 + omega * ddomain_ij;
				double domain_ij1 = domain[i][j];
				delta[i][j] = abs(domain_ij1 - domain_ij0);
			}
		}
		maxDelta = maxChange(delta);
		vec.push_back(maxDelta);
		//setBoundaryT(domain);
		setBoundaryP(domain);
		it++;
	}
	printMat2File("1020.csv", domain);
	printMatrix(domain, imax, jmax, "Point Gauss-Seidel with Over-Relaxation:");
	cout << "\nDelta = " << maxDelta << endl;
	cout << "\nIterations: " << it << endl;
	//printVector(vec);
	//printVec2File("DeltaVector_PGS_OR_2020_w15.csv", vec);
}

int main()
{
	double T[imax][jmax] = { {0} };			// Temperatue Matrix
	double P[imax][jmax] = { {0} };			// Pressure Matrix
	double S[imax][jmax] = { {0} };			// Source term
	double ExactT[imax][jmax] = { {0} };	// Exact Solution matrix
	double err[imax][jmax] = { {0} };		// Error matrix
	vector<double> delVec;					// Delta vector to store max change per iteration
	delVec.reserve(10000);					// Hopefully never more than 10000 iterations
	double delta[imax - 2][jmax - 2];		// Change matrix
	double tol = pow(10, -10);				// Change tolerance
	double omega = 1.5;						// Over-relaxation factor (w = 1 = No over-relaxation)
	double area = dx * dy;

	//exactSol(ExactT);
	//PGS_OR(T, delta, tol, delVec, omega, S);
	//error(err, T, ExactT);

	//printMat2File("error_2020_w15.csv", err);
	//printMatrix(ExactT, imax, jmax, "Exact");
	//printMat2File("ExactSolution_1010.csv", ExactT);
	//printMatrix(err, imax, jmax, "Error");

	source(S);
	PGS_OR(P, delta, tol, delVec, omega, S);

	double P0505 = 1.0 / 4 * (P[imax / 2][jmax / 2] + P[imax / 2 + 1][jmax / 2] + P[imax / 2][jmax / 2 + 1] + P[imax / 2 + 1][jmax / 2 + 1]);

	//cout << "\nL2 Norm = " << L2Norm(err) << endl;

	cout << "\nP(1/2,1/2) = " << P0505 << endl;

	//printMatrix(S, imax, jmax, "Source");
	
	getchar();
	return 0;
} 
