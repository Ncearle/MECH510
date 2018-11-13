%% MECH 510 - Assignment 3
% Nicholas Earle

clear; clc; close all

phi = 0:2*pi/100:2*pi;

%% Second Order Upwind Scheme - 1st Derivative

LU2_1 = @(phi) -(3 -4*cos(phi) + cos(2*phi) + 1i*(4*sin(phi) - sin(2*phi)));

LU2 = LU2_1(phi);

plot(1/2*real(LU2), 1/2*imag(LU2));
xlabel('Real \lambda *(u/\Deltax)');
ylabel('Imag \lambda *(u/\Deltax)');
title('Eigenvalues for 2nd order upwind and centred schemes');
grid on;

%% Second Order Centred Scheme - 2nd Derivative

LC2_2 = @(phi) cos(phi) - 1;

LC2 = LC2_2(phi);

hold on;
plot(2*real(LC2), 2*imag(LC2));


%% Combined

LUC2 = @(phi) -17 + 22*cos(phi) - 5*cos(2*phi) - 5i*(4*sin(phi) - sin(2*phi));

LUC = LUC2(phi);

plot(1/10*real(LUC), 1/10*imag(LUC));
legend('2nd Upwind - 1st Deriv.','2nd Centred - 2nd Deriv.', 'Combined');

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
g3 = 1 + z + 1/2*z.^2 + 1/6*z.^3;
g4 = 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4;

b = 1;
% Calculate magnitude of g
gmag2 = abs(g2);
gmag3 = abs(g3);
gmag4 = abs(g4);
% Plot contours of gmag
figure();
hold on;
contour(phi,y,gmag2,[b b], 'r','LineWidth', 2);
contour(phi,y,gmag3,[b b], 'g','LineWidth', 2);
contour(phi,y,gmag4,[b b], 'b','LineWidth', 2);
xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('Stability boundary for RK2, RK3, RK4');
legend('RK2','RK3','RK4');
axis tight;
grid on;

b = 0.2:0.2:1;
figure()
contour(phi,y,gmag4,b,'ShowText', 'on');
xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('RK4 Contour for |\sigma|');
axis tight;
grid on;

%% Max time step for RK4
figure()
contour(phi, y, gmag4, [1 1], '--','LineWidth', 2);
hold on;
plot(real(LUC/10),imag(LUC/10));
plot(0.75*real(LUC/10), 0.75*imag(LUC/10));
plot(0.5*real(LUC/10), 0.5*imag(LUC/10));
plot(0.633*real(LUC/10), 0.633*imag(LUC/10));
xlabel('Real \lambda\Deltat *(u/\Deltax)');
ylabel('Imag \lambda\Deltat *(u/\Deltax)');
title('Maximum stable time step for combined scheme with RK4');
legend('RK4 Stability','\Deltat = 1','\Deltat = 0.75','\Deltat = 0.5','\Deltat = 0.633');
grid on;


%% Max amp factor for funky 3 stage time advance

phi = pi/2:pi/200:pi;
f = @(phi) 1/2*(-3 + 4*cos(phi) - cos(2*phi) - 1i*(4*sin(phi) - sin(2*phi)));

sig = @(CFL) max(abs(1 + 97/75*CFL*f(phi) + 24/25*CFL^2*f(phi).^2 + 16/75*CFL^2*f(phi).^3));

cmin = fminsearch(sig, 0.3)
sig_max = sig(cmin)


