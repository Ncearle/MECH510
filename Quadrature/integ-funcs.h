#ifndef INTEG_FUNCS_H_
#define INTEG_FUNCS_H_

double exponential(const double x);
double bessel(const double x);
double sinsin(const double x, const double y);
double sinsinh(const double x, const double y);
void richardson(const double fine, const double medium, const double coarse, double* order, double* extrapolated);
double trapezoidal(double(*fun)(const double), const double lowerLimit, const double upperLimit, const int numIntervals);
double gauss(double(*fun)(const double), const double lowerLimit, const double upperLimit, const int numIntervals, const int numQuadInt);


#endif
