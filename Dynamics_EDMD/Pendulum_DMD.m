clear all; close all; clc
%This is just a function Test Bed
dt = 0.01;
duration = 30;
thresh = 0; % 0 is no truncation, if < 0 then sets singular value threshold, positive intergers determins truncation rank if possible
% Careful with set truncation, minimal testing was performed

%% Set structures for EDMDc

x0 = [.25;.25];
xf = [0;0];
mpcset.lqrQmod = 1; % lqr Q weight to adjust R/Q ratio
mpcset.lqrR =10;  

mpcset.R = 0; % Previous Control penalty -not used
mpcset.N = 20; % prediction horizonpy
mpcset.Nu = 20; % control horizon
mpcset.Q = [1;.1];  % 
% mpcset.Q = 10; % scalar Q
mpcset.Ru = 10; %Ctrl cost
pset.nvar = 2;  %initial # of states
pset.order = 2; % order of poly nomial fit
pset.sine = 0; % fourier order %% NOTE should be 2^n


%% Start Main
name = DMDPendulumtest(x0,xf,dt,duration,thresh,mpcset,pset);

%% Plot
% name =  plotend(name,duration,0);
% path2data = '../Data';
% load([path2data,'/',name,'/',[name,'Data.mat']])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Function

function [name] = DMDPendulumtest(x0,xf,dt,duration,thresh,mpcset,pset)
%Global Constants
g = -9.81;   %Gravity
mn = 1;       %Pendulum Mass
M = 5;       %Cart Mass
L = 2;       %Arm Length
d = 0;       %Damping (Currently Unused)
b = 1;       %
Ru = mpcset.Ru;  
tspan = 0.0:dt:duration;
x0(1) = x0(1)*pi;
x0(2) = x0(2)*pi;
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

%% File Management
var1 = num2str(x0(1));
var2 = num2str(x0(2));
var3 = num2str(round(mpcset.lqrR/mpcset.lqrQmod,2,'significant')); %round to 2 significant figures for the lqr Q/R ratio
var4 = num2str(Ru);    
time = clock;
now = [num2str(time(5)),'-',num2str(time(4)),'-',num2str(time(3)),'-',num2str(time(2))];
ModelName = 'Pendulum_DMD_Uncontrol';
ModelName1 = [ModelName,var1,'pi_',var2,'pi_',var3,'_',var4,'Ru_',now];

path2data = ['../Data/',ModelName1,'/']; mkdir(path2data)
path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)

%% Input System Function
f = @(t,x,u,p)pendulumeom(x,u);

% LQR Controller
A = [0 1; -L/g 0];
B = [0;1/(mn*L^2)];
Q = diag(mpcset.lqrQmod * mpcset.Q);
R = mpcset.lqrR;
gain = lqr(A,B,Q,R);
u=@(x)-gain*(x - xf);
[~,y1] = ode45(@(t,x)f(t,x,u(x),0),tspan,x0);
uvals = zeros(1,length(y1));
for k=1:length(y1) %capture lqr control val
     uvals(1,k) = - gain*(y1(k,:)'-xf);
end
[~,y2] = ode45(@(t,x)f(t,x,0,0),tspan,x0);
%% Run DMD
y1 = y1';
y2 = y2';

for i = 1:length(y1(1,:))
%     y1c(1,i) = wrapTo2Pi(y1(1,i));
    y1c(1,i) = y1(1,i);
    y1c(2,i) = y1(2,i);
%     y2c(1,i) = wrapTo2Pi(y2(1,i));
    y2c(1,i) = y2(1,i);
    y2c(2,i) = y2(2,i);
end

[Phi ,omega ,Atilde ,b,Xdmd] = myDMD(y2c,2,dt);
[Phic ,omegac ,Atildec ,bc,Xdmdc] = myDMD(y1c,2,dt);

Xdmd = real(Xdmd);
Xdmdc = real(Xdmdc);

% for i = 1:length(Xdmd(1,:)) 
%     Xdmd(1,i) = wrapTo2Pi(Xdmd(1,i));
%     Xdmdc(1,i) = wrapTo2Pi(Xdmdc(1,i));
% end
for i = 1:length(Xdmd(1,:)) 
    for k = 1:length(Xdmd(:,1))
        error1(k,i) = (y2c(k,i)-Xdmd(k,i))^2/y2c(k,i)^2;
        error2(k,i) = (y1c(k,i)-Xdmdc(k,i))^2/y1c(k,i)^2;
    end
end


%% Plotting

figure(1);
subplot(3,1,1)
plot(y2(1,1:end-1),'LineWidth',5,'Color','b')
hold on
plot(Xdmd(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')

subplot(3,1,2)
plot(y2(2,1:end-1),'LineWidth',5,'Color','b')
hold on
plot(Xdmd(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')

subplot(3,1,3)
plot(error1(1,:),'LineWidth',2,'Color','b')
hold on
plot(error1(2,:),'LineWidth',2,'Color','r')
hold off 
ylim([-.1 5])
legend('\theta error','\theta'' error')
xlabel('t - centiseconds')
ylabel('Error')


figure(2);
subplot(3,1,1)
plot(y1(1,1:end-1),'LineWidth',5,'Color','b')
hold on
plot(Xdmdc(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')

subplot(3,1,2)
plot(y1(2,1:end-1),'LineWidth',5,'Color','b')
hold on
plot(Xdmdc(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')

subplot(3,1,3)
plot(error2(1,:),'LineWidth',2,'Color','b')
hold on
plot(error2(2,:),'LineWidth',2,'Color','r')
hold off 
ylim([-.1 5])
legend('\theta error','\theta'' error')
xlabel('t - centiseconds')
ylabel('Error')

%% SAVE RESULTS

DataStore.y1 = y1;
DataStore.y2 = y2;
DataStore.A1 = Atilde;
DataStore.A2 = Atildec;
DataStore.Xdmd = Xdmd;
DataStore.Xdmd2 = Xdmdc;
DataStore.tspan = dt;
save([path2data,[ModelName1,'Data.mat']])
disp(['File saved at ../Data/',ModelName1])
saveas(figure(1),[path2data, 'States_no_control.fig'])
saveas(figure(1),[path2figs, 'States_no_control.png'])
saveas(figure(2),[path2data, 'States_w_control.fig'])
saveas(figure(2),[path2figs, 'States_w_control.png'])
name = ModelName1;
end

