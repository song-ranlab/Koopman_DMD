clear all; close all; clc
%This is just a function Test Bed
dt = 0.01;
duration = 30;
thresh = 0; % 0 is no truncation, if < 0 then sets singular value threshold, positive intergers determins truncation rank if possible
% Careful with set truncation, minimal testing was performed

%% Set structures for EDMDc

x0 = [-1;0;1.01;0.0];
xf = [1;0;pi;0];
mpcset.lqrQmod = 100; % lqr Q weight to adjust R/Q ratio
mpcset.lqrR = 1;  

mpcset.R = 0; % Previous Control penalty -not used
mpcset.N = 20; % prediction horizon
mpcset.Nu = 20; % control horizon
mpcset.Q = [20;1;10;1];  % 
% mpcset.Q = 10; % scalar Q
mpcset.Ru = 10; %Ctrl cost
pset.nvar = 4;  %initial # of states
pset.order = 2; % order of poly nomial fit
pset.sine = 0; % fourier order %% NOTE should be 2^n


%% Start Main
name = DMDInvPendulumtest(x0,xf,dt,duration,thresh,mpcset,pset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Function

function [name] = DMDInvPendulumtest(x0,xf,dt,duration,thresh,mpcset,pset)
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
ModelName = 'Pendulum_Inverted_EDMD_';
ModelName1 = [ModelName, var1, '_', var2,'_',var3,'pi_',var4,'pi_',var5,'_',var6,'Ru_',now];

path2data = ['../Data/',ModelName1,'/']; mkdir(path2data)
path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)

%% Input System Function
f = @(t,x,u,p)cartpend(x,mn,M,L,g,d,u);
A = [0 1 0 0;
    0 -d/M b* mn*g/M 0;
    0 0 0 1;
    0 -b*d/(M*L) -b*(mn+M)*g/(M*L) 0];
B = [0; 1/M; 0; b*1/(M*L)];
Q = diag(mpcset.lqrQmod * mpcset.Q);
R = mpcset.lqrR;
gain = lqr(A,B,Q,R);
u=@(x)-gain*(x - xf);
[~,z1] = ode45(@(t,x)f(t,x,u(x),0),tspan,x0);
uvals = zeros(1,length(z1));
for k=1:length(z1) %capture lqr control val
     uvals(1,k) = - gain*(z1(k,:)'-xf);
end
[~,z2] = ode45(@(t,x)f(t,x,0,0),tspan,x0);
%% Run DMD

y1 = poolData(z1,pset.nvar,3,0)'; 
y2 = poolData(z1,pset.nvar,0,4)'; 

% y1 = fft(z1');
[Phi ,omega ,Atilde ,b,Xdmd] = myDMD(y1,size(y1(:,1),1),dt);
[Phi2 ,omega2 ,Atilde2 ,bc,Xdmd2] = myDMD(y2,size(y2(:,1),1),dt);

y1 = real(y1);
Xdmd = real(Xdmd);
Xdmd2 = real(Xdmd2);

% for i = 1:length(Xdmd(1,:)) 
%     Xdmd(1,i) = wrapTo2Pi(Xdmd(1,i));
%     Xdmdc(1,i) = wrapTo2Pi(Xdmdc(1,i));
% end
for i = 1:length(Xdmd(1,:)) 
    for k = 1:length(Xdmd(:,1))
        error1(k,i) = abs(y2(k,i)-Xdmd(k,i))/y2(k,i);
        error2(k,i) = abs(y1(k,i)-Xdmd2(k,i))/y1(k,i);
    end
end


%% Plotting

figure(1);
subplot(5,1,1)
plot(y1(1,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X')
title('State 1')
title('EDMD 2nd Poly Order - State 1')

subplot(5,1,2)
plot(y1(2,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X''')
title('EDMD 2nd Poly Order - State 2')

subplot(5,1,3)
plot(y1(3,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(3,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')
title('EDMD 2nd Poly Order - State 3')

subplot(5,1,4)
plot(y1(4,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd(4,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')
title('EDMD 2nd Poly Order - State 4')

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
title('EDMD 2nd Poly Order - Error')


figure(2);
subplot(5,1,1)
plot(y2(1,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd2(1,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X')
title('EDMD 2nd Fourier Order - State 1')

subplot(5,1,2)
plot(y2(2,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd2(2,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('X''')
title('EDMD 2nd Fourier Order - State 2')

subplot(5,1,3)
plot(y2(3,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd2(3,1:end-1),'LineWidth',2,'Color','r')
hold off
% xlim([tspan(1) tspan(end)])
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta')
title('EDMD 2nd Fourier Order - State 3')

subplot(5,1,4)
plot(y2(4,1:end-1),'LineWidth',3,'Color','b')
hold on
plot(Xdmd2(4,1:end-1),'LineWidth',2,'Color','r')
hold off
legend('Original','DMD Reconstructed')
xlabel('t - centiseconds')
ylabel('\theta''')
title('EDMD 2nd Fourier Order - State 4')

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
title('EDMD 2nd Fourier Order - Error')

%% SAVE RESULTS

DataStore.y1 = y1;
DataStore.y2 = y2;
DataStore.A1 = Atilde;
DataStore.A2 = Atilde2;
DataStore.Xdmd = Xdmd;
DataStore.Xdmd2 = Xdmd2;
DataStore.tspan = dt;
save([path2data,[ModelName1,'Data.mat']])
disp(['File saved at ../Data/',ModelName1])
saveas(figure(1),[path2data, 'States_no_control.fig'])
saveas(figure(1),[path2figs, 'States_no_control.png'])
saveas(figure(2),[path2data, 'States_w_control.fig'])
saveas(figure(2),[path2figs, 'States_w_control.png'])
name = ModelName1;
end

