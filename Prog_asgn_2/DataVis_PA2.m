%% MECH 510 - Programming Assignment 2
% Nicholas Earle
% Data Visualisation

clear;clc;close all;

%% Exact Solution

xmax = 1;
tmax = 1;
n = 101;

T = @(x,t) sin(2*pi*(2*t-x));

[x, t] = meshgrid(linspace(0,xmax, n),linspace(0,tmax, n));

T_exact = T(x, t);

figure();
imagesc(flipud(T_exact));
colorbar;
set(gca, 'YDir', 'normal');

figure();
surf(x, t, T_exact);
