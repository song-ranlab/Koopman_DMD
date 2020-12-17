function [name,error,erroravg] = DMDpendEDMD(x0,xf,duration,item)
%% Script Prep
close all

%Global Constants
g = -10;   %Gravity
m = 1;       %Pendulum Mass
M = 5;       %Cart Mass
L = 2;       %Arm Length
d = 1;       %Damping (Currently Unused)
b = 1;       %Not 100% sure what this is but I assume its a linearized trig

thresh = 1e-4; 
noisemag = 0.0001;
%Predicted time, still needs to be fixed for DMDc
predict = 5; 
%Time Parameters
dt = 0.01;
% duration = 20;
tspan = 0.0:dt:duration;
%Misc Labeling
ver = ['v1p',num2str(item)];
var1 = num2str(x0(1));
var2 = num2str(x0(2));

var3 = num2str(duration);
x0(1) = x0(1)*pi;
x0(2) = x0(2)*pi;

%% File Management
time = clock;
now = [num2str(time(5)),'-',num2str(time(4)),'-',num2str(time(3)),'-',num2str(time(2))];
ModelName = 'Pendulum_EDMDc_EKF_';
ModelName1 = [ModelName, var1, 'pi_', var2,'pi_',var3,now];
path2data = ['../Data/',ModelName1,'/']; mkdir(path2data)
% path2figs = ['../Data/',ModelName1,'/']; mkdir(path2figs)

%% Other Parameters
%Input System Function
f = @(t,x,u,p)pendulumeom(x,u);
%Input (simple) Control Matrix
%B = [0; 1/M; 0; b*1/(M*L)];
%ODE
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
%LQR
Q = [1 0;
    0 1];
R = .0001;
%DMDc Modes


%% Unforced System
[t,y0] = ode45(@(t,x)f(t,x,0,0),tspan,x0,ode_options);
for k=1:length(y0)

    uvals(1,k) = 0;
end

%% Add Noise
 y = y0';
for i = 1: length(y0)
    y(1,i) = wrapTo2Pi(y(1,i) + noisemag * randn);
    y(2,i) = wrapTo2Pi(y(2,i) + noisemag * randn);
end


%% LQR Controller
% A = [0 1 0 0;
% 0 -d/M -m*g/M 0;
% 0 0 0 1;
% 0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0];
% 
% gain = lqr(A,B,Q,R);
% u=@(x)-gain*(x - xf);
% %[~,y1] = ode45(@(t,x)f(x,u(x)),tspan,x0);
% [~,y1] = ode45(@(t,x)f(t,x,u(x),0),tspan,x0);
% uvals = zeros(1,length(y1));
% for k=1:length(y1)
% %     uvals(1,k) = - gain*(y1(k,:)'-xf);
%     uvals(1,k) = 0
% end
% 
% 
% [n,blah] = size(y0');
y = y0';

%% EDMD 
[A, Xdmd,Y] = eDMDpoly(y, dt, thresh); %performs polynomial fit with auto
z = y(1,:);
obs = poolData(z,1,4,0);


%%
%kalman Filter
xkm1 = Xdmd(:,1);
Q = ;

for i = 1: length{tspan)
    
   %zbefore(k+1) = zest;
    H = [1:0];
    zest = H * xkm1 +0;
    %Iteration
    Pb = A * P * A' + Q;     %P(k+1)
    J =  H * Pb * H' + 1;    %holder variable R = 1
    K = Pb * H' * inv(J);    %kalman gainzzzz
    P = Pb - K * H * Pb;     %Post measurement P
    xk = A*xkm1 + K*(zest-H*A*xkm1);
    
end
%%
% %Kaiser's Run MPC
% %MPC Parameters
%  x_REF = xf;
% % x_REF = xliftf;
% options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
%     'MaxIterations',100);
% N=5;                %prediction horizon
% Nu = N;             %Control horizon
% mH=3;               %moving horizon
% Q3 = 5*eye(4);      %performance weights - size dependent on number of original states 
% R = 0;              %du weights
% Ru = 1;             %u weights
% QH = 5;             %Phi  weights
% p_EDMDc.sys = sysmodel_EDMDc;
% p_EDMDc.nvar = size(A2,1);
% p_EDMDc.polyorder = 4;
% p_EDMDc.usesine = 0; 
% p_EDMDc.Hfun = H;
% tic
% % [y2, u2, t2, r2] = runMPC(f,duration,dt,N,Nu,x0,@ObjectiveFCN_EDMDc,[],QH,R,Ru,options,x_REF,p_EDMDc);
% [y2, u2, t2, r2] = runMPC(f,duration,dt,N,Nu,x0,@ObjectiveFCN_EDMDc,[],QH,R,Ru,options,x_REF,p_EDMDc,Xdmdc);
% toc

%% LQR Attempt
%In synopsis it appears LQR cannot handle the size of the A matricies that
%EDMDcPoly is outputting so we are now moving to alternate linear control
%methods such as MPC
% eigKgain = lqr(A2,B2,Q2,R2);
% %------------------------------------------------------
%  lsys = sys((A2-B2*Kgain),B2,C2,D2,dt);
%  lsim(lsys,uvals,tspan);
% %------------------------------------------------------

% %[~,y2] = ode45(@(t,x)((A*x)-Kgain*(x-xliftf)), tspan, augx(:,1));
% [~,y2] = ode45(@(t,x)((A2-B2*Kgain)*(x-xliftf)), tspan, augx(:,1));
% uvals2 = zeros(1,length(y1));
% for k=1:length(y1)
%     uvals2(1,k) = - gain*(y1(k,:)'-xliftf);   %- our goal
% end
 %% SAVE RESULTS
 DataStore.y0 = y0;
%  DataStore.y1 = y1;
%  DataStore.y2 = y2;
 DataStore.u1 = uvals;
%  DataStore.u2 = u2;
%  DataStore.t2 = t2;
%  DataStore.r2 = r2;
 DataStore.tspan = dt;
 DataStore.Augx = augx;
 DataStore.Atilda = A2;
%  DataStore.Btilda = B2;
 %DataStore.XDMD = Xdmd;
 DataStore.XDMDC = Xdmdc;
%DataStore.Xp = Xp;
 save([path2data,[ModelName1,'Data.mat']])
%% Reconstruct using found U
% [~,ynew] = ode45(@(t,x)f(x,u2),tspan,x0);


%% Plot Results
% Prep for error plot
for i = 1:length(tspan)
    %y0check(i,:) = abs(y0((i+1),:)-y0k(:,i)');
    y1check(i,:) = abs(y1((i+1),:)-ynew(:,i)');
end

%% Controlled results
datit = sprintf('Inverted pendulum DMDc check for %d seconds',duration);
figtit = [datit,' Initial conditions:', var1, ' ', var2,' ',var3,'pi ',var4,'pi ',var5,'s ',ver];
Xerrorplot(y1',ynew',y1check',figtit);


%% Save Figure

% Save Snapshots
% F = getframe(gcf);
hmp = sprintf('DMD_Pendulum_Errorw%d.png',1);
saveas(figure(1),[path2data, hmp])
%%
name = ModelName1;
error = y1check(end,:);
erroravg = [mean(y1check(:,1)) mean(y1check(:,2)) mean(y1check(:,3)) mean(y1check(:,4))];
end

