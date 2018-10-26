clear;clc;

imax = 12;
jmax = 12;

T = zeros(jmax,imax);

for i = 1:imax
    T(1,i) = - T(2,i);
    T(jmax,i) = -T(jmax-1,i)+2*cos(pi*(i-1)/10);    
end
for j = 1:jmax
    T(j,1) = T(j,2);
    T(j,imax) = T(j,imax-1);
end


for j = 2:jmax-1
    for i = 2:imax-1
        T(i,j) = 1/4*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
    end
end

exact = zeros(10,10);
for j = 1:10
    for i = 1:10
        exact(j,i) = (cos(pi*(i-.5)/10)*sinh(pi*(j-.5)/10))/sinh(pi);
    end
end

        

        