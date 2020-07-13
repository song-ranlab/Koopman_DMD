function [name] = DMDquadrotor(x0,xf,duration)
%% Script Prep
close all

%Global Constants
g = 9.81;
m = .028;  %Vehicle Mass
b = 0.01; %base width
l=.025;%Arm Length
I_x = 2.3951 *10^-5; %from paper
I_y = 2.3951 *10^-5;
I_z = 1.8580 *10^-5;

Stat = [I_x; I_y; I_z; m];
%Predicted time, still needs to be fixed for DMDc
predict = 5; 
%Time Parameters
dt = 0.001;
% duration = 20;
tspan = 0.0:dt:duration;
%Misc Labeling
ver = 'path';
var1 = num2str(x0(1));
var2 = num2str(x0(2));
var3 = num2str(x0(3));
var4 = num2str(x0(7));
var5 = num2str(x0(8));
var6 = num2str(x0(9));
var7 = num2str(duration);
x0(7) = x0(7)*pi;
x0(8) = x0(8)*pi;
x0(9) = x0(9)*pi;
x0(10) = x0(10)*pi;
x0(11) = x0(11)*pi;
x0(12) = x0(12)*pi;
%% File Management
ModelName = 'Quadcopter_Trigfun';
ModelName1 = [ModelName, var1, '_', var2,'_',var3,'_',var4,'pi_',var5,'pi_',var6,'pi_',var7,'_',ver];
path2data = ['../Data/',ModelName1,'/']; mkdir(path2data)
% path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)

%% Other Parameters
%Input System Function
f = @(x,u)quadtrig(x,u,Stat);
%Input (simple) Control Matrix
% B = [0;0;0;(1/m);0;(1/m);0;0;0;(1/Ix);(1/Iy);(1/Iz)];
%ODE
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
%LQR
Q = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 10 0 0 0 0 0 0 0 0;
     0 0 0 0 10 0 0 0 0 0 0 0;
     0 0 0 0 0 10 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 0 0 10 0 0;
     0 0 0 0 0 0 0 0 0 0 10 0;
     0 0 0 0 0 0 0 0 0 0 0 1];
     
R = .01;
%DMDc Modes
r = 0;
p = 0;

%% Unforced System
[t,y0] = ode45(@(t,x)f(x,[0;0;0;0],Stat),tspan,x0,ode_options);
%% LQR Controller
A = [   0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 -g 0 0 0 0;
        0 0 0 0 0 0 g 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0];

B = [   0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 0 0 0; 
        0 0 0 0;
        1/m 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 I_x 0 0;
        0 0 I_y 0;
        0 0 0 I_z];


gain = lqr(A,B,Q,R);
[~,y1] = ode45(@(t,x)f(x,(-gain*(x-xf)),Stat),tspan,x0);
uvals = zeros(1,length(y1));
for k=1:length(y1)
    uvals(1,k) = - gain*(y1(k,:)'-xf);
end

%% DMD Prediction setup
% %Set learning time 
% lspan = ((duration-predict)/dt); 
% 
% for i = 1: lspan
%     y0l(i,1) = y0(i,1);
%     y0l(i,2) = y0(i,2);
%     y0l(i,3) = y0(i,3);
%     y0l(i,4) = y0(i,4);
%     y1l(i,1) = y1(i,1);
%     y1l(i,2) = y1(i,2);
%     y1l(i,3) = y1(i,3);
%     y1l(i,4) = y1(i,4);
% end

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
    y1check(i,:) = abs(y1((i+1),:)-Xp(:,i)');
end


%% Controlled results
datit = sprintf('Inverted pendulum DMDc check for %d seconds',duration);
figtit = [datit,' Initial conditions:', var1, ' ', var2,' ',var3,'pi ',var4,'pi ',var5,'s ',ver];

Xerrorplot(y1',Xp,y1check',figtit);
 %% SAVE RESULTS

 DataStore.y0 = y0;
 DataStore.y1 = y1;
 DataStore.u1 = uvals;
 DataStore.tspan = dt;
 DataStore.Atilda = Atild;
 DataStore.Btilda = Btild;
 DataStore.XDMD = Xdmd;
 DataStore.XDMDC = Xdmd2;
 DataStore.Xp = Xp;

%% Save Figure
save([path2data,[ModelName1,'Data.mat']])
% Save Snapshots
% F = getframe(gcf);
hmp = sprintf('DMD_Quadrotor_Errorw%d.png',1);
saveas(figure(1),[path2data, hmp])
%%
name = ModelName1;
error = y1check(end,:);
[x,y] = size(y1check');
for i = 1:y
    erroravg(i) = mean(y1check(:,i));
end
end

