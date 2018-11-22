%% MECH 510 - Assignment 3
% Nicholas Earle

clear; clc; close all

phi = 0:2*pi/100:2*pi;
Re = 50;    % Reynolds number
Pr = 0.7;   % Prandtl number
dx = 0.2;   % length of x cells
dy = 0.1;   % height of y cells
u = 3;      % value of ubar
v = 0;      % value of v

%% Second order centred scheme for energy equation

Lamx = @(phi) -u/dx * 1i*sin(phi) - 2/(Re*Pr*dx^2) * (1-cos(phi));
Lamy = @(phi) -v/dy * 1i*sin(phi) - 2/(Re*Pr*dy^2) * (1-cos(phi));

Lx = Lamx(phi);
Ly = Lamy(phi);

plot(real(Lx), imag(Lx));
hold on;
plot(real(Ly), imag(Ly));
xlabel('Real \lambda');
ylabel('Imag \lambda');
title('Eigenvalues for 2nd order centred scheme applied to energy equation');
legend('\lambda for X','\lambda for Y');
grid on;

%% Runge-Kutta Schemes

% Specify x range and number of points
x0 = -3;
x1 = 3;
Nx = 901;
% Specify y range and number of points
y0 = -3;
y1 = 3;
Ny = 901;
% Construct mesh
xv = linspace(x0,x1,Nx);
yv = linspace(y0,y1,Ny);
[phi,y] = meshgrid(xv,yv);

% Calculate z
z = phi + 1i*y;

g2 = 1 + z + 1/2*z.^2;

b = 1;
% Calculate magnitude of g
gmag2 = abs(g2);

% % Plot contours of gmag
% figure();
% hold on;
% contour(phi,y,gmag2,[b b], 'r','LineWidth', 2);
% xlabel('Real \lambda\Deltat');
% ylabel('Imag \lambda\Deltat');
% title('Stability boundary for RK2');
% legend('RK2');
% axis tight;
% grid on;

%% Max time step for RK2 in X
figure()
contour(phi, y, gmag2, [1 1], 'k--','LineWidth', 1);
hold on;

plot(0.07*real(Lx), 0.07*imag(Lx));
plot(0.065*real(Lx), 0.065*imag(Lx));
plot(0.0625*real(Lx), 0.0625*imag(Lx)); % Max Time Step "X"
plot(0.06*real(Lx), 0.06*imag(Lx));
plot(0.05*real(Lx), 0.05*imag(Lx));

xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('Maximum stable time step for energy equation with RK2 in x');
legend('RK2 Stability','\Deltat = 0.07','\Deltat = 0.065','\Deltat = 0.0625','\Deltat = 0.06','\Deltat = 0.05');
grid on;

%% Max time step for RK2 in Y
figure()
contour(phi, y, gmag2, [1 1], 'k--','LineWidth', 1);
hold on;

plot(0.2*real(Ly), 0.2*imag(Ly));
plot(0.175*real(Ly), 0.175*imag(Ly));   % Max Time step "Y"
plot(0.15*real(Ly), 0.15*imag(Ly));

xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('Maximum stable time step for energy equation with RK2 in y');
legend('RK2 Stability','\Deltat = 0.2','\Deltat = 0.175','\Deltat = 0.15');
grid on;