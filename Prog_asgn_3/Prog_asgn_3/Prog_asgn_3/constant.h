//==================================================
// Header file containing constants and includes
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <ctime>

using namespace std;

// Domain Constants
extern constexpr int imax = 102;
extern constexpr int jmax = 42;

extern constexpr double xmax = 5;
extern constexpr double ymax = 1;
extern constexpr double dx = xmax / (imax - 2);
extern constexpr double dy = ymax / (jmax - 2);

extern constexpr int ubar = 3;		// average velocity in x [m/s]
// extern constexpr double CFL = 0.001;	// CFL 80% of maximum for RK2
// extern constexpr double dt = CFL * dx / ubar;	//Time step

extern constexpr double dt = 0.004;

// Constants
extern constexpr double pi = M_PI;	// Pi
extern constexpr int Re = 50;		// Reynolds number
extern constexpr double Pr = 0.7;	// Prandtl number
extern constexpr double Ec = 0.1;	// Eckert number

// For Parts 1,2
extern constexpr int u0 = 1;		// Initial x velocity
extern constexpr int v0 = 1;		// Initial y velocity
extern constexpr int T0 = 1;		// Initial temperature