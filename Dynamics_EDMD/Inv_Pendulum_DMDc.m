clear all; close all; clc
%This is just a function Test Bed
dt = 0.001;
duration = 100;
thresh = 0; % 0 is no truncation, if < 0 then sets singular value threshold, positive intergers determins truncation rank if possible
% Careful with set truncation, minimal testing was performed

%% Set structures for EDMDc

x0 = [-1;0;1.01;0.0];
xf = [1;0;pi;0];
mpcset.lqrQmod = 1; % lqr Q weight to adjust R/Q ratio
mpcset.lqrR = 1000;  

mpcset.R = 0; % Previous Control penalty -not used
mpcset.N = 20; % prediction horizon
mpcset.Nu = 20; % control horizon
mpcset.Q = [20;1;10;1];  % 
% mpcset.Q = 10; % scalar Q
mpcset.Ru = 50; %Ctrl cost
pset.nvar = 4;  %initial # of states
pset.order = 2; % order of poly nomial fit
pset.sine = 0; % fourier order %% NOTE should be 2^n


%% Start Main
name = DMDcartpendEDMDPOLYtest(x0,xf,dt,duration,thresh,mpcset,pset);

%% Main Function

function [name] = DMDcartpendEDMDPOLYtest(x0,xf,dt,duration,thresh,mpcset,pset)
%Global Constants
g = -9.81;   %Gravity
mn = 1;       %Pendulum Mass
M = 5;       %Cart Mass
L = 2;       %Arm Length
d = 0;       %Damping (Currently Unused)
b = 1;       %
Ru = mpcset.Ru;  
tspan = 0.0:dt:duration;
x0(3) = x0(3)*pi;
x0(4) = x0(4)*pi;
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

%% File Management
var1 = num2str(x0(1));
var2 = num2str(x0(2));
var3 = num2str(x0(3));
var4 = num2str(x0(4));
var5 = num2str(round(mpcset.lqrR/mpcset.lqrQmod,2,'significant')); %round to 2 significant figures for the lqr Q/R ratio
var6 = num2str(Ru);    
time = clock;
now = [num2str(time(5)),'-',num2str(time(4)),'-',num2str(time(3)),'-',num2str(time(2))];
ModelName = 'Pendulum_Inverted_EDMDc_';
ModelName1 = [ModelName, var1, '_', var2,'_',var3,'pi_',var4,'pi_',var5,'_',var6,'Ru_',now];

path2data = ['../Data/',ModelName1,'/']; mkdir(path2data)
path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)

%% Input System Function
f = @(t,x,u,p)cartpend(x,mn,M,L,g,d,u);

% LQR Controller
A = [0 1 0 0;
    0 -d/M b* mn*g/M 0;
    0 0 0 1;
    0 -b*d/(M*L) -b*(mn+M)*g/(M*L) 0];
B = [0; 1/M; 0; b*1/(M*L)];
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
u1 = uvals;
%[Atilde, Btilde, nstates, augx] = eDMDCpoly(y, u1, dt, thresh,pset);
[Phi ,omega ,lambda ,modes, Xdmdc, Atilde, Btilde, r] = DMDcfin(y1, uvals, dt, thresh);
[Phi2 ,omega2 ,lambda2 ,modes2, Xdmd, Atilde2, Btilde2, r2] = DMDcfin(y2, zeros(size(uvals)), dt, 4);
Xdmd= real(Xdmd);

ynew(:,1) = x0;
for i = 2: length(tspan-1)
    ynew(:,i) = Atilde * ynew(:,i-1) + Btilde * uvals(i-1);
end
for i = 1:length(Xdmd(1,:)) 
    for k = 1:length(Xdmd(:,1))
        error1(k,i) = abs(y1(k,i)-ynew(k,i))/y1(k,i);
        error2(k,i) = abs(y2(k,i)-Xdmd(k,i))/y2(k,i);
    end
end


%% Plotting

figure(1);
subplot(5,1,1)
plot(y1(1,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(ynew(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X')
title('DMDc with Control - State 1')
subplot(5,1,2)
plot(y1(2,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(ynew(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X''')
title('DMDc with Control - State 2')
subplot(5,1,3)
plot(y1(3,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(ynew(3,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')
title('DMDc with Control - State 3')
subplot(5,1,4)
plot(y1(4,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(ynew(4,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')
title('DMDc with Control - State 4')
subplot(5,1,5)
plot(error1(1,:),'LineWidth',2,'Color','b')
hold on
plot(error1(2,:),'LineWidth',2,'Color','r')
hold on
plot(error1(3,:),'LineWidth',2,'Color','g')
hold on
plot(error1(4,:),'LineWidth',2,'Color','m')
hold off 
ylim([-.1 5])
legend('X Error','X'' Error','\theta error','\theta'' error')
xlabel('t - centiseconds')
ylabel('Error')
title('DMDc with Control - Error')

figure(2);
subplot(5,1,1)
plot(y2(1,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X')
title('DMDc with no Control - State 1')
subplot(5,1,2)
plot(y2(2,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X''')
title('DMDc with no Control - State 2')
subplot(5,1,3)
plot(y2(3,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(3,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')
title('DMDc with no Control - State 3')
subplot(5,1,4)
plot(y2(4,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(4,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')
title('DMDc with no Control - State 4')
subplot(5,1,5)
plot(error2(1,:),'LineWidth',2,'Color','b')
hold on
plot(error2(2,:),'LineWidth',2,'Color','r')
hold on
plot(error2(3,:),'LineWidth',2,'Color','g')
hold on
plot(error2(4,:),'LineWidth',2,'Color','m')
hold off 
ylim([-.1 5])
legend('X Error','X'' Error','\theta error','\theta'' error')
xlabel('t - centiseconds')
ylabel('Error')
title('DMDc with no Control - Error')
%% SAVE RESULTS

DataStore.y1 = y1;
DataStore.y2 = y2;
DataStore.u1 = u1;
DataStore.tspan = dt;
DataStore.Atilda = Atilde;
DataStore.Atilda2 = Atilde2;
save([path2data,[ModelName1,'Data.mat']])
disp(['File saved at ../Data/',ModelName1])
saveas(figure(1),[path2data, 'States_no_control.fig'])
saveas(figure(1),[path2figs, 'States_no_control.png'])
saveas(figure(2),[path2data, 'States_w_control.fig'])
saveas(figure(2),[path2figs, 'States_w_control.png'])
name = ModelName1;
end

%% Lifting for EDMDc
function [A, B, nstates, Y] = eDMDCpoly(StateData, InputData, dt, thresh,pset)
%Takes a full state, input vector and a time step and performs a polynomial
%regression and then DMDc on that new state. 
%Lift
gamma = [StateData; InputData];
Y = poolData(gamma',pset.nvar+1,pset.order,pset.sine)';       %Polynomial Fitting function from Kutz et Al

%% EDMD Reconstruct
[Phi ,omega ,lambda ,b, Xdmdc,A,B,r] = DMDcfin(Y, InputData, dt, thresh);
nstates = size(A,1);
end
