//#define _USE_MATH_DEFINES
//#include <iostream>
//#include <math.h>
//#include <fstream>
//#include <vector>
//
//using namespace std;
//
//constexpr int imax = 11;
//constexpr double xmax = 1.0;
//constexpr double dx = xmax / (imax - 1);
//
//constexpr double u = 2.0;
//constexpr double CFL = 0.4;
//constexpr double dt = CFL * dx / u;
//constexpr double tend = 0.1;
//constexpr int num_t_steps = tend / dt;
//constexpr int tmax = num_t_steps*2;			// Double rows for the intermediate step in RK2
//
//
//void printArray(double arr[imax]) // Prints a vector to the terminal
//{
//	for (int i = 0; i < imax; i++)
//	{
//		cout << arr[i] << "\t";
//	}
//	cout << endl;
//}
//
//// Initial Conditions
//void setInit(double *arr)
//{
//	for (int i = 0; i < imax; i++)
//	{
//		arr[i] = -sin(2 * M_PI*(i - 0.5) / (imax - 1));
//	}
//}
//
//// Boundary Conditions
////void setBoundary(double(&domain)[tmax][imax], int t_step) 
////{
////		domain[t_step][0] = 2 * sin(4 * M_PI * i * dt/2) - domain[t_step][1];
////}
//
//// 2nd order upwind flux integral
//double *FI2U(double *arr) 
//{
//	double temp[imax];
//	for (int i = 2; i < imax; i++) 
//	{
//		temp[i] = (3 * arr[i] - 4 * arr[i - 1] + arr[i - 2]) / (2 * dx);
//	}
//	return temp;
//}
//
//// Runge Kutta 2 time advance
////void RK2(vector<double> &vec, void(*scheme)(vector<double>))		
////{
////	setInit(vec);
////	printVector(vec);
////
////	vector<double> w_int(imax, 0.0);
////
////	for (int i = 0; i < tmax; i++)
////	{
////		// Intermediate step
////		w_int = vec + dt/2 * 
////
////	}
////
////
////
////}
//
//// 1st order upwind flux integral
//vector<double> FI1U(vector<double> &vec) 
//{
//	vector<double> temp(imax, 0.0);
//	for (int i = 1; i < vec.size(); i++) 
//	{
//		temp[i] = (vec[i] - vec[i-1]) / dx;
//	}
//	return temp;
//}
//
//// Explicit Euler time advance
////void EE(vector<vector<double> > &vec) 
////{
////	vector<double> temp(imax, 0);
////	setInit(temp);
////
////	for (int t_step = 1; t_step < tmax; t_step++) 
////	{
////		vector<double> temp;
////		temp = FI1U(vec);
////	}
////}
//
//
//int main()
//{
//
//	double T[imax] = { 0 };
//	setInit(T);
//	double *U = FI2U(T);
//	double wint[imax];
//	for (int i = 0; i < imax; i++)
//	{
//		double temp[imax];
//		temp[i] = dt / 2 * U[i];
//		wint[i] = T[i] + temp[i];
//	}
//
//	printArray(T);
//	cout << endl;
//	printArray(U);
//	cout << endl;
//	printArray(wint);
//
//	getchar();
//	return 0;
//}