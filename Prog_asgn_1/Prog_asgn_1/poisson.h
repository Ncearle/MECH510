#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

template<typename T>
void printMatrix(T mat, const int N, const int M, const char* S);
void setBoundary(double(&domain));
double exactSol(double(&exact));
double error(double(&error), double(&domain), double(&exact));
