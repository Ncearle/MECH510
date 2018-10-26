%% MECH 510 - Programming Assignment 1
% Data visualization
% Nicholas Earle

clear;clc;

%% Part 1

del_2020_w10 = fscanf(fopen('DeltaVector_PGS_OR_2020_w10.csv', 'r'), '%f');
del_2020_w13 = fscanf(fopen('DeltaVector_PGS_OR_2020_w13.csv', 'r'), '%f');
del_2020_w15 = fscanf(fopen('DeltaVector_PGS_OR_2020_w15.csv', 'r'), '%f');

T_1010_w10 = csvread('Data_PGS_OR_1010_w10.csv');
T_1010_w15 = csvread('Data_PGS_OR_1010_w15.csv');
T_2020_w10 = csvread('Data_PGS_OR_2020_w10.csv');
T_2020_w13 = csvread('Data_PGS_OR_2020_w13.csv');
T_2020_w15 = csvread('Data_PGS_OR_2020_w15.csv');
T_4040_w15 = csvread('Data_PGS_OR_4040_w15.csv');
T_8080_w15 = csvread('Data_PGS_OR_8080_w15.csv');

exact_1010 = csvread('ExactSolution_1010.csv');
exact_2020 = csvread('ExactSolution_2020.csv');
exact_4040 = csvread('ExactSolution_4040.csv');
exact_8080 = csvread('ExactSolution_8080.csv');

error_1010_w10 = csvread('error_1010_w10.csv');
error_1010_w15 = csvread('error_1010_w15.csv');
error_2020_w10 = csvread('error_2020_w10.csv');
error_2020_w13 = csvread('error_2020_w13.csv');
error_2020_w15 = csvread('error_2020_w15.csv');
error_4040_w15 = csvread('error_4040_w15.csv');
error_8080_w15 = csvread('error_8080_w15.csv');

fclose('all');

figure(1);
plot(1:length(del_2020_w10),del_2020_w10, 1:length(del_2020_w13), del_2020_w13, 1:length(del_2020_w15), del_2020_w15);
title('Max Change vs Iteration for 20x20 matrix with w = 1');
xlabel('Iteration');
ylabel('Max Change');
legend('w = 1.0','w = 1.3','w = 1.5');

%Iterations for tol = 1e-7

its_1010_w10 = 210;
its_1010_w15 = 93;

its_2020_w10 = 637;
its_2020_w13 = 410;
its_2020_w15 = 290;

its_4040_w10 = 1888;
its_4040_w15 = 876;

its_8080_w10 = 5854;
its_8080_w15 = 2548;

% Tolerance set to 1e-10
its_1010_w15_10 = 142;
its_2020_w15_10 = 487;
its_4040_w15_10 = 1649;
its_8080_w15_10 = 5546;


%Tolerance set to 1e-7
L2Norm_1010_w10 = 0.0023427;
L2Norm_2020_w10 = 0.000623219;
L2Norm_4040_w10 = 0.000160659;
L2Norm_8080_w10 = 5.91651e-05;


% L2 Norm and iterations for different meshes, w = 1.5, tol = 1e-10
L2Norm_1010_w15_10 = 0.00234271;
L2Norm_2020_w15_10 = 0.000623228;
L2Norm_4040_w15_10 = 0.000159244;
L2Norm_8080_w15_10 = 4.01541e-05;


figure(2);
plot([10,20,40,80],[L2Norm_1010_w15_10, L2Norm_2020_w15_10, L2Norm_4040_w15_10, L2Norm_8080_w15_10]);

figure();
imagesc(flipud(exact_1010)); 
colorbar;
set(gca, 'YDir', 'normal');

figure();
imagesc(flipud(exact_2020));
colorbar;
set(gca, 'YDir', 'normal');

figure();
imagesc(flipud(exact_4040));
colorbar;
set(gca, 'YDir', 'normal');

figure();
imagesc(flipud(exact_8080));
colorbar;
set(gca, 'YDir', 'normal');

%% Part 2 - Application

P_1010 = csvread('P_1010.csv');
P_2020 = csvread('P_2020.csv');
P_4040 = csvread('P_4040.csv');
P_8080 = csvread('P_8080.csv');

Phalf_1010 = 4.81406;
Phalf_2020 = 4.88991;
Phalf_4040 = 4.91645;
Phalf_8080 = 4.92758;

Pits_1010 = 391;
Pits_2020 = 1311;
Pits_4040 = 4590;
Pits_8080 = 16470;

% Error Calculation
phi1 = Phalf_8080;
phi2 = Phalf_4040;
phi3 = Phalf_2020;

h1 = sqrt(1/6400);
h2 = sqrt(1/1600);
h3 = sqrt(1/400);

r21 = h2/h1;
r32 = h3/h2;

ep21 = phi2-phi1;
ep32 = phi3-phi2;

p = 1/log(r21)*abs(log(abs(ep32/ep21)));

phi21_ext = (r21*phi1-phi2)*(r21-1);
phi32_ext = (r32*phi3-phi2)*(r32-1);

e21_a = abs((phi1-phi2)/phi1);

e21_ext = abs((phi21_ext-phi1)/phi21_ext);

GCI21_fine = 1.25*e21_a/(r21-1);

figure();
imagesc(flipud(P_1010)); 
colorbar;
set(gca, 'YDir', 'normal');
title('Pressure 10x10 mesh');

figure();
imagesc(flipud(P_2020)); 
colorbar;
set(gca, 'YDir', 'normal');
title('Pressure 20x20 mesh');

figure();
imagesc(flipud(P_4040)); 
colorbar;
set(gca, 'YDir', 'normal');
title('Pressure 40x40 mesh');

figure();
imagesc(flipud(P_8080)); 
colorbar;
set(gca, 'YDir', 'normal');
title('Pressure 80x80 mesh');

