#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

using namespace std;

template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S);
void printVector(vector<double> vec);

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