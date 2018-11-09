#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>

double bisection(double(*f)(double), double lowLim, double upLim, double tol) {

	double err = INFINITY;
	double a = lowLim, b = upLim;
	double c, root;
	double it = 0;

	while (err > tol) {
		
		it++;
		std::cout << "Iteration: " << it << "\t\ta = " << a << "\t\tb = " << b << "\n";
		
		c = (a + b) / 2;
		double fa = f(a), fb = f(b), fc = f(c);

		if (fa * fc < 0) {
			b = c;
		}
		else {
			a = c;
		}

		err = abs(a - b);
	}

	if (f(a) == 0) {
		root = a;
	}
	else if (f(b) == 0) {
		root = b;
	}
	else {
		root = (a + b) / 2;
	}

	std::cout << "\n\nA root is: " << root;
	std::cout << "\n\nIterations = " << it<<"\n\n";

	return root;
}

double secant(double(*f)(double), double lowLim, double upLim, double tol) {
	
	double err = INFINITY;
	double a = lowLim, b = upLim;
	double c, m, root;
	double it = 0;

	while (err > tol) {
		it++;
		std::cout << "Iteration: " << it << "\t\ta = " << a << "\t\tb = " << b << "\n";
		
		m = (f(b) - f(a)) / (b - a);
		c = b - f(b) / m;

		a = b;
		b = c;

		err = abs(a - b);
	}

	if (f(a) == 0) {
		root = a;
	}
	else if (f(b) == 0) {
		root = b;
	}
	else {
		root = (a + b) / 2;
	}

	std::cout << "\n\nA root is: " << root;
	std::cout << "\n\nIterations = " << it << "\n\n";

	return root;

}

double newton(double(*f)(double), double(*df)(double), double x0, double tol) {

	double delta = INFINITY;
	double it = 0;

	while (delta > tol) {
		double xnew = x0 - f(x0) / df(x0);
		delta = abs(xnew - x0);
		x0 = xnew;
		std::cout << x0 << "\n";
		it++;
	}

	double root = x0;

	std::cout << "\n\nA root is: " << root;
	std::cout << "\n\nIterations = " << it << "\n\n";

	return root;
	
}
double testFun1(double x) {
	return (x - 1)*(x - 2)*(x - 3);
}
double dxFun1(double x) {
	return (x - 2)*(x - 3) + (x - 1)*(x - 3) + (x - 1)*(x - 2);
}

double testFun2(double x) {
	return (pow(x, -4) - 1);
}
double dxFun2(double x) {
	return -4 * pow(x, -5);
}

double testFun3(double x) {
	return 20 * (x - M_PI)*exp(-25 * (pow(x - M_PI, 2)));
}
double dxFun3(double x) {
	return 20 * (1 + (x - M_PI)*(-50 * (x - M_PI))) * exp(-25 * pow(x - M_PI, 2));
}

int main() {

	//bisection(testFun1, 0, M_PI, pow(10, -5));
	//bisection(testFun2, 0.5, M_PI, pow(10, -5));

	//secant(testFun1, 0, M_PI, pow(10, -5));
	//secant(testFun2, 0.5, M_PI, pow(10, -5));

	//newton(testFun1, dxFun1, 4, pow(10, -5));
	//newton(testFun2, dxFun2, 4, pow(10, -5));
	newton(testFun3, dxFun3, 0, pow(10, -5));

	getchar();
	return 0;
}