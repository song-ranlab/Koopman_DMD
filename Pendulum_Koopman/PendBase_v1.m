%% Script Prep
clear all, close all, clc
%Input System Function
f = @(x,u)SimpPend(x,u);
%Input (simple) Control Matrix
B = [0;1];
%Initial Conditions
thetai = 0.5;
thetadi = 0.5;
x0 = [thetai;thetadi];
%Terminal Conditions
xf = [pi;0];
%Predicted time, still needs to be fixed for DMDc
predict = 5; 
%Time Parameters
dt = 0.01;
duration = 20;
tspan = 0.0:dt:duration;
%Misc Labeling
ver = 'v1p1';
var1 = num2str(x0(1));
var2 = num2str(x0(2));
var3 = num2str(duration);
%% File Management
ModelName = 'Pendulum_Simple_';
ModelName1 = [ModelName, var1, '_', var2,'_',var3,'_',ver];
path2data = ['../Data/',ModelName1]; mkdir(path2data)
path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)
%% Other Parameters
%ODE
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
%LQR
Q = [ 10 0; 0 10];
R = 1;
%DMDc Modes
r = 2;
p = 2;

%% Unforced System
[t,y0] = ode45(@(t,x)f(x,0),tspan,x0,ode_options);
%% LQR Controller
A = [0 1; -1 0];

gain = lqr(A,B,Q,R);
[~,y1] = ode45(@(t,x)f(x,(-gain*(x-xf))),tspan,x0);
uvals = zeros(1,length(y1));
for k=1:length(y1)
    uvals(1,k) = - gain*(y1(k,:)'-xf);
end

%% DMD Prediction setup
%Set learning time 
lspan = ((duration-predict)/dt); 

for i = 1: lspan
    y0l(i,1) = y0(i,1);
    y0l(i,2) = y0(i,2);
end

%% DMD 
[Mode,ceval,deval,magmode,Xdmd] = DMD(y0',r,dt); 
y0k = real(Xdmd);

%% DMDc 
[Mode2,ceval2,deval2,magmode2,Xdmd2,Atild,Btild,Xp] = DMDcii(y1',uvals,r,p,dt);
y1k = real(Xdmd2);

%% Plot Results
% Prep for error plot
for i = 1:length(y0k)
    y0check(i,:) = abs(y0((i+1),:)-y0k(:,i)');
    y1check(i,:) = abs(y1((i+1),:)-y1k(:,i)');
end
datit = sprintf('Simple pendulum DMD check for %d seconds',duration);
figure('NumberTitle', 'off', 'Name', datit);
%title('Inverted pendulum DMD check for %d seconds',duration)
subplot(2,1,1)
plot(y0check(:,1))
hold on 
plot(y0k(1,2:length(y0k))','x')
hold on
plot(y0(:,1))
hold off
title('Uncontrol Error - Theta')
legend ('Error','DMD','Normal')
ylabel('\theta');
xlabel('t(ms)');

subplot(2,1,2)
plot(y0check(:,2))
hold on 
plot(y0k(2,2:length(y0k)'),'x')
hold on
plot(y0(:,2))
hold off
title('Uncontrol Error - Theta_{dot}')
legend ('Error','DMD','Normal')
ylabel('\theta''');
xlabel('t(ms)');

datit = sprintf('Simple pendulum DMDc check for %d seconds',duration);
figure('NumberTitle', 'off', 'Name', datit);
%title('Inverted pendulum DMD check for %d seconds',duration)
subplot(2,1,1)
plot(y1check(:,1))
% plot(y1(:,1))
hold on 
plot(y1k(1,2:length(y1k))','x')
% plot(y1(:,1))
hold on
plot(y1(:,1))
hold off
title('Control Error - Theta')
legend ('Error','DMD','Normal')
ylabel('\theta');
xlabel('t(ms)');

subplot(2,1,2)
plot(y1check(:,2))
hold on 
plot(y1k(2,2:length(y1k)'),'x')
hold on
plot(y1(:,2))
hold off
title('Uncontrol Error - Theta_{dot}')
legend ('Error','DMD','Normal')
ylabel('\theta''');
xlabel('t(ms)');
 
 %% SAVE RESULTS

 DataStore.y0 = y0;
 DataStore.y1 = y1;
 DataStore.u1 = uvals;
 DataStore.tspan1 = dt;
 DataStore.Atilda = Atild;
 DataStore.Btilda = Btild;
 DataStore.XDMD = Xdmd;
 DataStore.XDMDC = Xdmd2;
% 
save([path2data,[ModelName1,'Data.mat']])
% Save Snapshots
 for i = 1:2
    hmp = sprintf('DMD_Pendulum_Error%d.png',i);
    saveas(figure(1),[path2figs, hmp])
 end
