% This code is made by Ahmed Eltahan

%* This code intends to solve 1st order ODE Runge–Kutta–Fehlberg procedure which is 6th order accuracy 


clc; close all; clear all;

%% Inputs and Declaration
h = 0.01; % solution stepsize
x = 0:h:pi/4; % input vector


% Function declaration
f = @(x, y) (y-x-1)^2+2;

% Initial conditions
y_RKF(1) = 1;


gamma = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]; % vector related to RKF procdure

%% Core: Runge Kutta 6th order procedure
for i = 1:length(x)-1

k1 = h*f(x(i), y_RKF(i));                                                       k1 = double(k1);
k2 = h*f(x(i)+h/4, y_RKF(i)+k1/4);                                              k2 = double(k2);
k3 = h*f(x(i)+3/8*h, y_RKF(i)+3/32*k1+9/32*k2);                                 k3 = double(k3);
k4 = h*f(x(i)+12/13*h, y_RKF(i)+1932/2197*k1-7200/2197*k2+7296/2197*k3);        k4 = double(k4);
k5 = h*f(x(i)+h, y_RKF(i)+439/216*k1-8*k2+3680/513*k3-845/4104*k4);             k5 = double(k5);
k6 = h*f(x(i)+h/2, y_RKF(i)-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5);   k6 = double(k6);

K = [k1, k2, k3, k4, k5, k6];

y_RKF(i+1) = y_RKF(i) + sum(K.*gamma); % new solution

end

%% Analytical Exact Solution
y_exact = tan(x) + x + 1;
error_RKF = y_exact - y_RKF; % error between 


%% Plot and Comparison
w = figure(1);
plot(x, y_exact, x, y_RKF, 'LineWidth', 2);
grid on
xlabel('x')
ylabel('y')
title('Plot of the Exact and Numerical RKF Solutions')
legend('RKF solution', 'Exact solution')
print(w, '-dpng', '-r720', 'Exact_VS_RKF')

w = figure(2);
plot(x, error_RKF, 'LineWidth', 2);
grid on
xlabel('x')
ylabel('Error RKF')
title('Plot of the Error between Exact and Numerical RKF Solutions')
print(w, '-dpng', '-r720', 'Error_RKF')
