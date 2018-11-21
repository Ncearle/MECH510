//==================================================
// Header file containing print functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "error_fcns.h"
#include "constant.h"

void printVec(vector<double> vec);		// Prints a vector to the console
void printVec2D(vector<vector<double>> vec);		// Prints a 2D vector to the console
void vec2D2File(const char* fileName, vector<vector<double>> &vec); // Prints a 2D vector to a file
void vec1D2File(const string fileName, vector<double> &vec); // Prints 1D vector to a file

void printVec(vector<double> vec) 
{
		for (int i = 0; i < vec.size(); i++)
		{
			cout << fixed << setprecision(4) << vec[i] << "  ";
		}
		cout << endl;
}

void printVec2D(vector<vector<double>> vec) 
{
	for (int j = vec.size()-1; j > -1; j--)
	{
		for (int i = 0; i < vec[j].size(); i++)
		{
			cout << fixed << setprecision(3) << vec[j][i] << "  ";
		}
		cout << endl;
	}
	cout << endl;
}


void vec2D2File(const string fileName, vector<vector<double>> &vec) // Prints 2D vector to a file
{			
	ofstream outfile;
	outfile.open(fileName);
	for (int j = vec.size() - 2; j > 0; j--)
	{
		for (int i = 1; i < vec[j].size() - 1; i++)
		{
			outfile << vec[j][i] << " ";		// Separated by a space for .dat files
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}

void vec1D2File(const string fileName, vector<double> &vec) // Prints 1D vector to a file
{
	ofstream outfile;
	outfile.open(fileName);
	for (int i = 0; i < vec.size(); i++)
	{
		outfile << vec[i] << " ";	// Separated by a space for .dat files
	}
	outfile << endl;
	outfile.close();
}