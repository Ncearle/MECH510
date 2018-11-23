%% Data Fitting

clear; clc; close all;

% [x,y] = importdata('dataset1.dat', 2);

x = [1.05, 1.259133, 1.459230, 1.695443, 1.873318, 2.055083, 2.326705, 2.503156, 2.742026, 2.868842];
y = [0.705794, 1.174020, 1.716892, 2.170531, 2.642141, 3.067350, 3.665453, 4.049727, 4.634206, 4.899706];
sig = 0.02;


Sxy = 0;
Sxx = 0;
Sx = 0;
S = 0;
Sy = 0;

N = length(x);

for n = 1:N
    Sxx = Sxx + x(n)^2/sig^2;
    Sx = Sx + x(n)/sig^2;
    Sxy = Sxy + x(n)*y(n)/sig^2;
    S = S + 1/sig^2;
    Sy = Sy + y(n)/sig^2;
end


delta = S*Sxx - Sx^2;
a1 = (S*Sxy - Sx*Sy)/delta;
a2 = (Sxx*Sy - Sx*Sxy)/delta;
yxn = a1 * x + a2;
plot(x, y, '*', x, yxn);
