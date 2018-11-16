%% MECH 510 - Programming Assignment 3 - Data visualization
% Nicholas Earle

clear; clc; close all;

% 10 x 10 domain
T_10 = csvread('T_10.csv');
u_10 = csvread('u_10.csv');
v_10 = csvread('v_10.csv');
ErrFI_10 = csvread('ErrFI_10.csv');
ErrS_10 = csvread('ErrS_10.csv');

L2_FI10 = 0.0956385;
L2_S10 = 0.000577231;

xv10 = linspace(0.05, 0.95, 10);
yv10 = xv10;
[x10, y10] = meshgrid(xv10, yv10);

% 20 x 20 domain
T_20 = csvread('T_20.csv');
u_20 = csvread('u_20.csv');
v_20 = csvread('v_20.csv');
ErrFI_20 = csvread('ErrFI_20.csv');
ErrS_20 = csvread('ErrS_20.csv');

L2_FI20 = 0.0337088;
L2_S20 = 0.000146137;

xv20 = linspace(0.05, 0.95, 20);
yv20 = xv20;
[x20, y20] = meshgrid(xv20, yv20);

% 40 x 40 domain
T_40 = csvread('T_40.csv');
u_40 = csvread('u_40.csv');
v_40 = csvread('v_40.csv');
ErrFI_40 = csvread('ErrFI_40.csv');
ErrS_40 = csvread('ErrS_40.csv');

L2_FI40 = 0.0213721;
L2_S40 = 3.66494e-05;

xv40 = linspace(0.05, 0.95, 40);
yv40 = xv40;
[x40, y40] = meshgrid(xv40, yv40);

figure();
surf(x10, y10, ErrFI_10);
figure();
surf(x20, y20, ErrFI_20);
figure();
surf(x40, y40, ErrFI_40);

