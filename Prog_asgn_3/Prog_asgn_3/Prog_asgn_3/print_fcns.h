//==================================================
// Header file containing print functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "error_fcns.h"
#include "constant.h"

void printVec2D(vector<vector<double>> vec);		// Prints a 2D vector to the console
void vec2File(const char* fileName, vector<vector<double>> &vec); // Prints a 2D vector to a file

void printVec2D(vector<vector<double>> vec) 
{
	for (int j = jmax - 1; j > -1; j--)
	{
		for (int i = 0; i < imax; i++)
		{
			cout << fixed << setprecision(4) << vec[j][i] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}


void vec2File(const char* fileName, vector<vector<double>> &vec) // Prints 2D vector to a file
{			
	ofstream outfile;
	outfile.open(fileName);
	for (int j = jmax - 2; j > 0; j--)
	{
		for (int i = 1; i < imax - 1; i++)
		{
			outfile << vec[j][i] << " ";		// Separated by a space for .dat files
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}