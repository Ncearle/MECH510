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
colorbar;
title('Flux Integral Error (80 x 80 mesh)');
xlabel('x'); ylabel('y'); zlabel('T');
figure();
surf(x80, y80, ErrS_80);
shading interp;
colorbar;
title('Source Error (80 x 80 mesh)');
xlabel('x'); ylabel('y'); zlabel('T');

%% Part 7

figure(3);
G1 = load('Grad_x100_w1.dat');
d1 = diff(G1);
l1 = find(d1 < 1e-4, 1);
x1 = l1/length(G1)*100;
plot(G1);
title('\DeltaT_{y0} for W = 1; X_{max} = 100');
xlabel('Cells x');
ylabel('\DeltaT');

figure(2);
G05 = load('Grad_x50_w0.dat');
d05 = diff(G05);
l05 = find(d05 < 1e-4/0.5, 1);
x05 = l05/length(G05)*50;
plot(G05);
title('\DeltaT_{y0} for W = 0.5; X_{max} = 50');
xlabel('Cells x');
ylabel('\DeltaT');

figure(1);
G025 = load('Grad_x25_w025.dat');
d025 = diff(G025);
l025 = find(d025 < 1e-4/0.25, 1);
x025 = l025/length(G025)*25;
plot(G025);
title('\DeltaT_{y0} for W = 0.25; X_{max} = 25');
xlabel('Cells x');
ylabel('\DeltaT');

figure(4);
G2 = load('Grad_x200_w2.dat');
d2 = diff(G2);
l2 = find(d2 < 1e-4/2, 1);
x2 = l2/length(G2)*200;
plot(G2);
title('\DeltaT_{y0} for W = 2; X_{max} = 200');
xlabel('Cells x');
ylabel('\DeltaT');

figure(5);
G4 = load('Grad_x600_w4.dat');
d4 = diff(G4);
l4 = find(d4 < 1e-4/4, 1);
x4 = l4/length(G4)*600;
plot(G4);
title('\DeltaT_{y0} for W = 4; X_{max} = 600');
xlabel('Cells x');
ylabel('\DeltaT');

figure(6);
G8 = load('Grad_x1200_w8.dat');
d8 = diff(G8);
l8 = find(d8 < 1e-4/8, 1);
x8 = l8/length(G8)*1200;
plot(G8);
title('\DeltaT_{y0} for W = 8; X_{max} = 1200');
xlabel('Cells x');
ylabel('\DeltaT');


