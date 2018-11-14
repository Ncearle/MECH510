%% Data Fitting

clear; clc; close all;

[x,y] = importdata('dataset1.dat', '\t', 2);




Sxy = 0;
Sxx = 0;
Sx = 0;
S = 0;
Sy = 0;
for n = 1:N
Sxx = Sxx + x(n)^2/sig(n)^2;
Sx = Sx + x(n)/sig(n)^2;
Sxy = Sxy + x(n)*y(n)/sig(n)^2;
S = S + 1/sig(n)^2;
Sy = Sy + y(n)/sig(n)^2;
end


delta = S*Sxx - Sx^2;
a1 = (S*Sxy - Sx*Sy)/delta;
a2 = (Sxx*Sy - Sx*Sxy)/delta;
yxn = a1 * x + a2;
plot(x, y, x, yxn);
