//==================================================
// Header file containing constants and includes
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

vector<vector<double>> exactFlux()
{
	vector<vector<double>> ExFI(jmax, vector<double>(imax));
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			ExFI[j][i] = u0 * T0*pi*cos(2 * pi*x)*y*sin(pi*y) + v0 * T0*pi*x*cos(pi*x)*cos(2 * pi*y) + 2 * T0*pow(pi, 2)*cos(pi*x)*sin(pi*y) / (Re*Pr);
		}
	}
	return ExFI;
}

vector<vector<double>> exactSource()
{
	vector<vector<double>> ExS(jmax, vector<double>(imax));
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			ExS[j][i] = Ec / Re * (2 * pow(u0*pi*cos(pi*x)*y, 2) + 2 * pow(v0*pi*x*sin(pi*y), 2) + pow(u0*sin(pi*x) + v0 * cos(pi*y), 2));
		}
	}
	return ExS;
}