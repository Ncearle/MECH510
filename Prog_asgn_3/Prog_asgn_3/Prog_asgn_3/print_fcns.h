//==================================================
// Header file containing all of the print functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include <iomanip>
#include "error_fcns.h"

using namespace std;

template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S);		// Prints a matrix to the console
void printVec2D(vector<vector<double>> vec);									// Prints a vector to the console
void printM2F(const char* fileName, double(&matrix)[jmax][imax]);
void printMat2File(const char* fileName, double(&matrix)[jmax][imax]);	// Prints a matrix to a file

template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S) {
	cout << "\n" << S << ":\n";
	for (int i = N - 2; i > 0; i--) 
	{
		for (int j = 1; j < N - 1; ++j)
		{
			cout << fixed << setprecision(6) << *(*(mat + i) + j) << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void printVec2D(vector<vector<double>> vec) 
{
	for (int j = jmax - 2; j > 0; j--)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			cout << fixed << setprecision(6) << vec[j][i] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}


void printM2F(const char* fileName, double(&matrix)[jmax][imax]) // Prints Matrix to a file
{			
	ofstream outfile;
	outfile.open(fileName);
	for (int j = jmax - 2; j > 0; j--)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			outfile << matrix[j][i] << ", ";
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}

void printMat2File(const char* fileName, double(&matrix)[jmax][imax]) // Prints Matrix to a file
{			
	ofstream outfile;
	outfile.open(fileName);
	for (int j = jmax - 2; j > 0; j--)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			outfile << matrix[j][i] << ", ";
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}