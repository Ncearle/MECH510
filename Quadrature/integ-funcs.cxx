#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "integ-funcs.h"

double exponential(const double x)
{
	return exp(-x*x/2);
}

double bessel(const double x)
{
	return _j0(x)*_y1(x);
}

double sinsin(const double x, const double y)
{
	return sin(M_PI*x)*sin(M_PI*y)*(M_PI*M_PI)/4;
}

double sinsinh(const double x, const double y)
{
	return sin(M_PI*x)*sinh(M_PI*y);
}

void richardson(const double fine, const double medium, const double coarse, double * order, double * extrapolated)
{
	double k = log2 ((coarse - medium) / (medium - fine));
	double ue = fine - (medium - fine) / (pow(2, k) - 1);
	std::cout << "\n\nRichardson:\n\n";
	std::cout << "k = " << k << "\nue = " << ue<<std::endl;
	*order = k;
	*extrapolated = ue;
}

double trapfun(const double x)
{
	return M_PI / 2 * sin(M_PI*x);
}

double trapezoidal(double(*fun)(const double), const double lowerLimit, const double upperLimit, const int numIntervals)
{
	std::cout << "\n\nTrapezoidal:\n\n";
	double sum = 0;
	double n = (upperLimit - lowerLimit) / numIntervals;

	for (int i = (lowerLimit + n)/n; i < (upperLimit - n)/n; i++) {
		sum += trapfun(i*n);
		/*std::cout<<"\n"<<sum;*/
	}

	double I = n / 2 * (trapfun(lowerLimit) + 2*sum + trapfun(upperLimit));
	std::cout << "\nI = " << I;
	return I;
}

double gauss(double(*fun)(const double), const double lowerLimit, const double upperLimit, const int numIntervals, const int numQuadInt)
{
	std::cout << "\n\nGauss Quadrature:\n\n";
	double n = (upperLimit - lowerLimit) / numIntervals;

	/*switch (numQuadInt) {
		case 1:

		case 2:

		case 3:
	}*/
	
	return 0;
}

int main()
{
	double k, ue;
	richardson(4.01541e-05, 1.5924e-04, 6.2323e-04, &k, &ue);

	/*double tc, tm, tf;
	tc = trapezoidal(trapfun, 0, 1, 10);
	tm = trapezoidal(trapfun, 0, 1, 20);
	tf = trapezoidal(trapfun, 0, 1, 40);

	richardson(tf, tm, tc, &k, &ue);*/

	getchar();
	return 0;
}
