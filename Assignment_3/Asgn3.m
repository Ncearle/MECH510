%% MECH 510 - Assignment 3
% Nicholas Earle

clear;clc;

phi = 0:2*pi/100:2*pi;

%% Second Order Upwind Scheme - 1st Derivative

LU2_1 = @(phi) -(3 -4*cos(phi) + cos(2*phi) + sqrt(-1)*(sin(phi) - sin(2*phi)));

LU2 = LU2_1(phi);

plot(real(LU2), imag(LU2));
xlabel('Real \lambda');
ylabel('Imag \lambda');
title('Eigenvalues for 3rd order upwind-biased and upwind schemes');
grid on;

%% Second Order Centred Scheme - 2nd Derivative

LC2_2 = @(phi) cos(phi) - 1;

LC2 = LC2_2(phi);

hold on;
plot(real(LC2), imag(LC2));
legend('Second Order Upwind Scheme - 1st Derivative','Second Order Centred Scheme - 2nd Derivative');

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

% figure()
% contour(phi,y,gmag3,[1 1],'ShowText', 'on');
% xlabel('Real \lambda\Deltat');
% ylabel('Imag \lambda\Deltat');
% title('Contour for |\sigma| = 1');
% grid on;
