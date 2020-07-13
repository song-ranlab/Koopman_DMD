function [y] = ODE4HB(F_xy, tspan, x0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = tspan(2)-tspan(1);
x = tspan;
y = zeros(4,length(tspan)); 
x0(3) = wrapTo2Pi(x0(3));
y(:,1) = x0;


for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y(:,i));
    k_2 = F_xy(x(i)+0.5*h,y(:,i)+0.5*h*k_1);
    k_3 = F_xy((x(i)+0.5*h),(y(:,i)+0.5*h*k_2));
    k_4 = F_xy((x(i)+h),(y(:,i)+k_3*h));

    y(:,i+1) =y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
     y(3,i+1) = wrapTo2Pi(y(3,i+1));
%     y(4,i+1) = wrapTo2Pi(y(4,i+1));
end
end

