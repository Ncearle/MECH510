%% MECH 510 - Programming Assignment 2
% Nicholas Earle
% Data Visualisation

clear; clc; close all;

%% Exact Solution

xmax = 1;
imax = 41;
dx = xmax/(imax-1);

tmax = 1;
u = 2.0;
CFL = 0.4;
dt = CFL * dx / u;
n = int8(tmax/dt);

%% Exact Solution
T = @(x,t) sin(2*pi*(2*t-x));

[x, t] = meshgrid(linspace(0, xmax, 1/dx), linspace(0, tmax, 1/dt+1));

T_exact = T(x, t);

figure();
imagesc(flipud(T_exact));
colorbar;
set(gca, 'YDir', 'normal');
title('Exact Solution with 40 cells')
xlabel('Cell');
ylabel('Timestep');

figure();
surf(x, t, T_exact);

%% Numerical Solution;

% x = 20, CFL = 0.4, 2nd order upwind
T_20_04 = csvread('Data_RK2_FI2U_20_04.csv');       % Temperature
E_20_04 = csvread('Error_RK2_FI2U_20_04.csv');      % Error
L2N_20_04 = csvread('L2Norm_RK2_FI2U_20_04.csv');   % L2 norm

T_40_04 = csvread('Data_RK2_FI2U_40_04.csv');       % Temperature
E_40_04 = csvread('Error_RK2_FI2U_40_04.csv');      % Error
L2N_40_04 = csvread('L2Norm_RK2_FI2U_40_04.csv');   % L2 norm

% for i = 1:n
%    figure(3);
%    plot(x(i,:), T_20_04(101-i,:), 'LineWidth', 2);
%    axis([0 1 -1 1]);
%    hold on;
%    plot(x(i,:), T_exact(101-i,:), 'LineWidth', 2);
%    hold off;
%    drawnow;
%    
%    figure(4);
%    plot(x(i,:), E_20_04(101-i,:), 'LineWidth', 2);
%    axis([0 1 -1 1]);
%    drawnow;
% end


x20 = linspace(0, xmax, 20/xmax);
x40 = linspace(0, xmax, 40/xmax);
t20 = linspace(0, tmax, 100);
t40 = linspace(0, tmax, 200);

figure();
plot(t20, L2N_20_04,t40, L2N_40_04);
title('L2 norm for each time step');
xlabel('Time');
ylabel('L2norm');
legend('20 cells','40 cells');

figure();
plot(x20, E_20_04(end,:), x40, E_40_04(end,:));
title('Error at t = 1');
xlabel('x');
ylabel('Error (num. - exact)');
legend('20 cells','40 cells');

%% Part 2 Stability

x500 = linspace(0, 20+1/500, 501);
T_500_100 = csvread('Data_RK2_FI2U_500_100.csv');
T_500_101 = csvread('Data_RK2_FI2U_500_101.csv');
T_500_102 = csvread('Data_RK2_FI2U_500_102.csv');

T_500_1007_30 = csvread('Data_RK2_FI2U_500_1007.csv');
T_500_1008_30 = csvread('Data_RK2_FI2U_500_1008.csv');
T_500_1009_30 = csvread('Data_RK2_FI2U_500_1009.csv');


figure()
plot(x500, T_500_102(end,:));
hold on
plot(x500, T_500_101(end, :))
plot(x500, T_500_100(end, :))
legend('\Deltat = 0.0102','\Deltat = 0.0101','\Deltat = 0.0100');
title('Stability of RK2 with 2nd order upwind scheme, for different timesteps');
xlabel('X');
ylabel('T');
axis([0 20 -2 2]);

figure()
plot(x500, T_500_1009(end,:));
hold on
plot(x500, T_500_1008(end, :))
plot(x500, T_500_1007(end, :))
legend('\Deltat = 0.01009','\Deltat = 0.01008','\Deltat = 0.01007');
title('Stability of RK2 with 2nd order upwind scheme, for different timesteps');
xlabel('X');
ylabel('T');
axis([0 20 -2 2]);

%% Part 3 More Schemes

Linf_RK2_FI2U_80 = csvread('LInf_RK2_FI2U_80.csv');
L1N_RK2_FI2U_80 = csvread('L1Norm_RK2_FI2U_80.csv');   

Linf_RK2_FI1U_80 = csvread('LInf_RK2_FI1U_80.csv');
L1N_RK2_FI1U_80 = csvread('L1Norm_RK2_FI1U_80.csv');   

Linf_EE_FI1U_80 = csvread('LInf_EE_FI1U_80.csv');
L1N_EE_FI1U_80 = csvread('L1Norm_EE_FI1U_80.csv');   

Linf_EE_FI2U_80 = csvread('LInf_EE_FI2U_80.csv');
L1N_EE_FI2U_80 = csvread('L1Norm_EE_FI2U_80.csv');   

t400 = linspace(0, 1, 400);
t200 = linspace(0, 1, 200);

figure();
plot(t400, Linf_RK2_FI2U_80, t200, Linf_RK2_FI1U_80, t400, Linf_EE_FI2U_80, t200, Linf_EE_FI1U_80);
legend('RK2 - 2U','RK2 - 1U','EE - 2U','EE - 1U');
title('Linf for each timestep for the different schemes');
xlabel('time');
ylabel('Linf');
figure();
plot(t400, L1N_RK2_FI2U_80, t200, L1N_RK2_FI1U_80, t400, L1N_EE_FI2U_80, t200, L1N_EE_FI1U_80);
legend('RK2 - 2U','RK2 - 1U','EE - 2U','EE - 1U');
title('L1 norm for each timestep for the different schemes');
xlabel('time');
ylabel('L1 norm');

figure();
plot(t400, Linf_RK2_FI2U_80, t400, Linf_EE_FI2U_80, t200, Linf_EE_FI1U_80);
legend('RK2 - 2U','EE - 2U','EE - 1U');
title('Linf for each timestep for the different schemes');
xlabel('time');
ylabel('Linf');
figure();
plot(t400, L1N_RK2_FI2U_80, t400, L1N_EE_FI2U_80, t200, L1N_EE_FI1U_80);
legend('RK2 - 2U','EE - 2U','EE - 1U');
title('L1 norm for each timestep for the different schemes');
xlabel('time');
ylabel('L1 norm');

%% L norms mesh refinement

L1_EE_1 = [0.0761122, 0.049345, 0.0279237, 0.0148212, 0.00763068, 0.00387097];
L2_EE_1 = [0.104468, 0.064482, 0.0357291, 0.0187928, 0.0096357, 0.004877859];
LI_EE_1 = [0.207714, 0.125343, 0.0688273, 0.0361063, 0.0184779, 0.00934554];

L1_EE_2 = [0.190377, 0.141569, 0.0849228, 0.045823, 0.0236799, 0.0120207];
L2_EE_2 = [0.208511, 0.161413, 0.0967912, 0.0520454, 0.0268401, 0.0136099];
LI_EE_2 = [0.298844, 0.259047, 0.157408, 0.0843862, 0.0433867, 0.0219611];

L1_RK2_1 = [0.324644, 0.216225, 0.129157, 0.0710064, 0.037308, 0.019134];
L2_RK2_1 = [0.398859, 0.268085, 0.160863, 0.0889079, 0.0468334, 0.0240463];
LI_RK2_1 = [0.726831, 0.498777, 0.304775, 0.169351, 0.0894109, 0.0459649];

L1_RK2_2 = [0.14749, 0.056309, 0.0165801, 0.00440982, 0.00113013, 0.000285575];
L2_RK2_2 = [0.188604, 0.0704658, 0.0207985, 0.00552953, 0.00141581, 0.000357553];
LI_RK2_2 = [0.328291, 0.155082, 0.0493493, 0.0133178, 0.0034302, 0.000864263];

cells = [10, 20, 40, 80, 160, 320];

figure();
loglog(cells, L1_RK2_2, cells, L1_RK2_1, cells, L1_EE_2, cells, L1_EE_1);
legend('RK2-2U','RK2-1U','EE-2U','EE-1U');
title('L1 norm for different mesh sizes');
xlabel('Number of cells');
ylabel('L1 norm');


figure();
loglog(cells, L2_RK2_2, cells, L2_RK2_1, cells, L2_EE_2, cells, L2_EE_1);
legend('RK2-2U','RK2-1U','EE-2U','EE-1U');
title('L2 norm for different mesh sizes');
xlabel('Number of cells');
ylabel('L2 norm');

figure();
loglog(cells, LI_RK2_2, cells, LI_RK2_1, cells, LI_EE_2, cells, LI_EE_1);
legend('RK2-2U','RK2-1U','EE-2U','EE-1U');
title('LInf norm for different mesh sizes');
xlabel('Number of cells');
ylabel('LInf norm');
