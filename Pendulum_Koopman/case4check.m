clear all, close all, clc

ModelName = 'Pendulum_';
ModelName1 = [ModelName, 'Uncontrol_m1_M5_0_0_p75pi_0_c4_'];
ModelName2 = [ModelName, 'Uncontrol_Hamiltonian_m1_M5_0_0_p75pi_0_c4_'];
path2data = ['../Data/',ModelName1]; mkdir(path2data)

path2figs = ['../Figures/PENDULUM/',ModelName1,'/']; mkdir(path2figs)
path2figs2 = ['../Figures/PENDULUM/',ModelName2,'/']; %mkdir(path2figs)

g = -9.81;
m = 1;  %Pendulum Mass
M = 5;  %Cart Mass
L = 2;%Arm Length
C = 0;  %Damping - Currently Unused
duration = 5;
dt = 0.001;
% Parameters
tspan = 0.0:dt:duration;
x0 = [0;0;0.000001;.00001]; % Initial Conditions

Q  = 1;     
R  = 1;    

%% Plot Settings
ivp=[x0(3); x0(4); g; M; L ; C];
arrow = 'ShowArrow';
interval = [0, duration];
% movie = true;
%%
    H = @(x)(.5* M * x(2)^2 +...
        .5 * m*((x(2)^2 + ...
        2*x(2)*L*x(4)*cos(x(3))+ ...
        (L*x(4)*sin(x(3)))^2)+...
        (L*x(4)*cos(x(3)))^2) -...
        m*g*L*(1-cos(x(3))));

gradH = @(x)([  0;...
                (M+m)*x(2)+m*L*x(4)*sin(x(3));
                m*L*x(2)*x(4)*cos(x(3))-m*g*L*sin(x(3));
                m*L*x(2)*sin(x(3))+m*L^2*x(4)]);

[U,V] = meshgrid([-2*pi:0.01:2*pi], [-2*pi:0.01:2*pi]);

[Y,X] = meshgrid([-2*pi:0.01:2*pi],[-4:0.01:4] );

%% Unforced
B=[0;1;0;0];
d=0;

f = @(t,x,u) ([x(2);
               ((1/(m*L*L*(M+m*(1-cos(x(3))^2))))*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + 0*m*L*L*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))*u + u);
           x(4);
           (1/(m*L*L*(M+m*(1-cos(x(3))^2))))*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - 0*m*L*cos(x(3))*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))]);
               
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

% % KRONIC: Case 4
H_ref = H([0,0,.1,0])%REF = -1;
H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-H_ref);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[y4] = ode4(@(t,x)f(t,x,-gain(x)*H_err(x)),tspan,x0);
%[y4] = ode4(@(t,x)f(t,x,0),tspan,x0);
uvals4 = zeros(1,length(y4)); for k=1:length(y4), uvals4(1,k) =-gain(y4(k,:))*H_err(y4(k,:)); end
[Hvals4,Jvals4] = evalCostFun_KoopEfun(H,y4,uvals4,Q,R,H_ref);

% Store results
DataStore.y4 = y4;
% DataStore.u4 = uvals4;
DataStore.H4 = Hvals4;
DataStore.J4 = Jvals4;
DataStore.tspan4 = tspan;

%% SAVE RESULTS
save([path2data,'/',[ModelName1,'Data.mat']])