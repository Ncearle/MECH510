%% MECH 510 - Programming Assignment 3 - Data visualization
% Nicholas Earle

clear; clc; close all;

%% Part 1, 2
% 10 x 10 domain
T_10 = importdata('T_10.dat');
ErrFI_10 = importdata('ErrFI_10.dat');
ErrS_10 = importdata('ErrS_10.dat');

L2_FI10 = 0.0877319;
L2_S10 = 0.000577231;

xv10 = linspace(0.05, 0.95, 10);
yv10 = xv10;
[x10, y10] = meshgrid(xv10, yv10);

% 20 x 20 domain
T_20 = importdata('T_20.dat');
ErrFI_20 = importdata('ErrFI_20.dat');
ErrS_20 = importdata('ErrS_20.dat');

L2_FI20 = 0.0223009;
L2_S20 = 0.000146137;

xv20 = linspace(0.05, 0.95, 20);
yv20 = xv20;
[x20, y20] = meshgrid(xv20, yv20);

% 40 x 40 domain
T_40 = importdata('T_40.dat');
ErrFI_40 = importdata('ErrFI_40.dat');
ErrS_40 = importdata('ErrS_40.dat');

L2_FI40 = 0.00559842;
L2_S40 = 3.66494e-05;

xv40 = linspace(0.05, 0.95, 40);
yv40 = xv40;
[x40, y40] = meshgrid(xv40, yv40);

% 80 x 80 domain
T_80 = importdata('T_80.dat');
ErrFI_80 = importdata('ErrFI_80.dat');
ErrS_80 = importdata('ErrS_80.dat');

L2_FI80 = 0.00140106;
L2_S80 = 9.16957e-06;

xv80 = linspace(0.05, 0.95, 80);
yv80 = xv80;
[x80, y80] = meshgrid(xv80, yv80);


figure();
surf(x10, y10, ErrFI_10);
shading interp;
figure();
surf(x20, y20, ErrFI_20);
shading interp;
figure();
surf(x40, y40, ErrFI_40);
shading interp;
figure();
surf(x80, y80, ErrFI_80);
shading interp;

%% Part 3-6



