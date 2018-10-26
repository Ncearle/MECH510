%% MECH 510 - Assignment 2
% Nicholas Earle

clear;clc;

x = 0:2*pi/100:2*pi;

%% Third Order Upwind-Biased Scheme

Lubk = @(phi) 3 + 4*cos(phi) + cos(2*phi) + sqrt(-1)*(-8*sin(phi) - sin(2*phi));

Lub = -Lubk(x);

plot(real(Lub), imag(Lub));
xlabel('Real \lambda');
ylabel('Imag \lambda');
title('Eigenvalues for 3rd order upwind-biased and upwind schemes');
grid on;

%% Third Order Upwind only Schee

Luok = @(phi) 11 - 18*cos(phi) + 9*cos(2*phi) - 2*cos(3*phi) + sqrt(-1)*(18*sin(phi) - 9*sin(2*phi) + 2*sin(3*phi));

Luo = Luok(x);

hold on;
plot(real(Luo), imag(Luo));
legend('3rd Order Upwind-biased','3rd Order Upwind');

%% Three-stage Runge-Kutta

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
[x,y] = meshgrid(xv,yv);
% Calculate z
z = x + 1i*y;
g3 = 1 + z + 1/2*z.^2 + 1/6*z.^3;

b = 0.1:0.1:1;
% Calculate magnitude of g
gmag = abs(g3);
% Plot contours of gmag
figure();
contour(x,y,gmag,b,'ShowText', 'on');
%axis([x0,x1,y0,y1]);
%axis('square');
xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('Contours for stability boundary for RK3');
legend('|\sigma|');
grid on;

figure()
contour(x,y,gmag,[1 1],'ShowText', 'on');
xlabel('Real \lambda\Deltat');
ylabel('Imag \lambda\Deltat');
title('Contour for |\sigma| = 1');
grid on;
