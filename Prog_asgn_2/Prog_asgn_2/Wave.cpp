#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include<vector>

using namespace std;

constexpr int imax = 11
constexpr double xmax = 1.0
constexpr double dx = xmax / (imax - 1)

constexpr double u = 2.0
constexpr double CFL = 0.4
constexpr double dt = CFL*dx/u;
constexpr double tend = 1;
constexpr int num_t_steps = tend/dt;


template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S) {
	cout << "\n" << S << ":\n";
	for (int i = N-2; i > 0; i--) {
		for (int j = 1; j < N-1; ++j)
			cout << *(*(mat + i) + j) << "  ";
		cout << endl;
	}
	cout << endl;
}

void printVector(vector<double> vec) {			// Prints a vector to the terminal
	for (int i = 0; i < vec.size(); i++)
	{
		cout << vec[i] << endl;
	}
}

// Initial Conditions
void setInitX(double(&domain)[num_t_steps][imax]) {
	for (int i = 1; i < imax; i++)
	{
		domain[0][i] = -sin(2*M_PI*(i+0.5)/(imax-1));
	}
}

// Boundary Conditions
void setBoundary(double(&domain)[num_t_steps][imax], int t_step){

	domain[t_step][0] = 2*sin(4*M_PI*tend/t_step);
}


int main(){

	double T[t_steps][imax] = {{0}};


	setInitX(T);

	printMatrix(T,num_t_steps, imax, "T" );

	getchar();
	return 0;
}
