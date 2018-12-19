function [c] = fit(RXPdata1)
D=length(RXPdata1);           %   Max meaurement distance
%
d1=(1:D);   %   1m to max measurment distance (D) in Increments (Inc)
% RXPdata1=[-45.0983; -53.3746; -54.4132; -56.8286; -59.2905; -60.2743;... 
%    -60.6919; -59.4938; -64.0525; -62.6163]; %    Measured data in dB
%
semilogy(d1,(RXPdata1(:,1)),'-or')          %    Plot measured data
hold on
%
%                %   Graph Set up
%
title ('Plot showing measured data','FontSize',14, 'FontWeight','bold')
xlabel('Distance (m) [Log scale]','FontSize',10, 'FontWeight','bold')
ylabel('Recieved Power (dB)','FontSize',10, 'FontWeight','bold')
%
hold on
grid on
%
              %   Regression
%                
logd1 = log10(d1');
A     = [ones(size(RXPdata1)),logd1];
c = A\RXPdata1
y = A*c;
plot(d1,y,'o-', d1,y,'*-')
%
legend ('Measured Data','Lin. Regression')