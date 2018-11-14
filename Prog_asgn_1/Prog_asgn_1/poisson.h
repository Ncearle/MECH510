#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

using namespace std;

template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S);
void printVector(vector<double> vec);
void printMat2File(const char* fileName, double(&matrix)[imax][jmax]);
void printVec2File(const char* fileName, vector<double> vec);


template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S) {
	cout << "\n" << S << ":\n";
	for (int i = N - 2; i > 0; i--) {
		for (int j = 1; j < N - 1; ++j)
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

void printMat2File(const char* fileName, double(&matrix)[imax][jmax]) {			// Prints Matrix to a file
	ofstream outfile;
	outfile.open(fileName);
	for (int i = imax - 2; i > 0; i--)
	{
		for (int j = 1; j < jmax - 1; j++)
		{
			outfile << matrix[i][j] << ", ";
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}

void printVec2File(const char* fileName, vector<double> vec) {			// Prints vector to a file
	ofstream outfile;
	outfile.open(fileName);
	for (int i = 0; i < vec.size(); i++)
	{
		outfile << vec[i] << endl;
	}
	outfile << endl;
	outfile.close();
}