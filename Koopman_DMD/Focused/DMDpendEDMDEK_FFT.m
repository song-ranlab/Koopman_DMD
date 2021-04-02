function [name,error,erroravg] = DMDpendEDMDEK(x0,xf,duration,item)
%% Script Prep
close all

%Global Constants
g = -10;   %Gravity
m = 1;       %Pendulum Mass
M = 5;       %Cart Mass
L = 2;       %Arm Length
d = 1;       %Damping (Currently Unused)
b = 1;       %Not 100% sure what this is but I assume its a linearized trig

thresh = 0; 
order = 16;
noisemag = 0.000;
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
name = ModelName1;

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
%% FFT Dis Stuff
% % Defineing Space
% L = duration * dt;          %Length of space
% N = size(tspan);            %# of points
% dx = L/N;
% x_fft = -L/2:dx:L/2-dx;     %# Y domain
% f_fft = 
% ut=fft2(y0);
% [t,utsol]=ode45(@(t,x)f_fft(t,x,0,0),tspan,fft(x0),ode_options);

% for j=1:size(B,1);  % Compute row-wise FFT
%     Cshift(j,:) = fftshift(fft(y0(j,:)));
%     C(j,:) = (fft(y0(j,:)));
% end
% 
% for j=1:size(C,2);  % Compute column-wise FFT
%     D(:,j) = fft(C(:,j));
% end



%% Add Noise
 y = y0';    %Change to proper SS form
for i = 1: length(y0)
    
    yn(2,i) = wrapToPi(y(2,i) + noisemag * randn);
    if i == 1
        yn(1,i) = y(1,i);
    else
        yn(1,i) = yn(2,i-1) * dt + yn(1,i-1);
    end
    yobs(1,i) = wrapToPi(y(1,i) + noisemag * randn);
end
% fftlift = fft2(yn);
%% EDMD 
%performs fourier fit with auto
[A_raw, Xdmd,Y] = eDMDFFT(yn, dt, order, thresh);
% z = y(1,:);
% zest = fft(yobs,order/2,1);% order/2 part is probably incorrect
zest = fft(yobs);
[r, o2] = size(Xdmd);

xlift = real(Xdmd);
Ylift = real(Y);
%if rank(A_raw) == rank(Y) 

    for i = 1:r
        for k = 1:o2
            errstate(i,k) = xlift(i,k) - Ylift(i,k);
        end
    end

    %% Check error in lifted state

    figure
    for i = 1:r
    %     if i == 10
    %         plot(errstate(i,:),'+')
    %     elseif i == 11
    %         plot(errstate(i,:),'v') 
    %     else
            plot(errstate(i,:))
    %     end
        hold on 
    end
    title('Error in lifted approximation')
    legend
%end

%% Check error in norm state

normregstate = real(ifft(Y));
normregstate(1,:) = wrapToPi(normregstate(1,:));
normtrunkstate = real(ifft(Xdmd));%% dont wrap theta dot
normtrunkstate(1,:) = wrapToPi(normtrunkstate(1,:));

% normregstate = wrapToPi(real(ifft(Y,order,1)));
% normtrunkstate = wrapToPi(real(ifft(Xdmd,order,1)));
% 
% figure
% for v = 1:2
%     check(v,:) = abs((normtrunkstate(v,1:length(tspan))-normregstate(v,1:length(tspan)-1))/normregstate(v,1:length(tspan)-1))*100;
%     plot(check(v,:))
%     hold on 
% end
% title('Error in returned state')
% legend('\theta','\theta''')
%%Plot States
t = 1:25:length(tspan)-1;
figure
% for z = 1:2
    plot(t,yn(1,t),'--r+')
    hold on
%     pause
    plot(t,yn(2,t),'--bx')
    hold on
%     pause
% end
t = 1:30:length(tspan)-1;
% for q =1:2
    plot(t,normregstate(1,t)','-.gv');
    hold on
%     pause
    plot(t,normregstate(2,t)','--mo');
    hold on
%     pause
% end
t = 1:15:length(tspan)-1;
% for v = 1:2
    plot(t,normtrunkstate(1,t)','-.r*');
    hold on 
%     pause
    plot(t,normtrunkstate(2,t)',':ks');
% end
title('Check States')
legend('Vanilla \theta','Vanilla \theta''','Lifted \theta','Lifted \theta''','Truncated \theta','Truncated \theta''')
%% kalman Filter

xkm1 = Ylift(:,1);
xk(:,1) = xkm1;
A = real(A_raw);
Q = .5*eye(size(A));

sig = sqrt(noisemag*.9); %for now need to find better value
R = sig^2;% * eye(size(A));
% Make H matrix
% H = zeros(order/2,order);
% for i = 1:order/2
%    H(i,i) = 1;
% end
H =[1 0];
% P = A*sig*A' + B*Q*B';
P = A*sig*A';
%Pcheck(1
 for i = 1:length(tspan)
    
   %zbefore(k+1) = zest;

%    zest = H * xkm1 +0;
    %Iteration
    Pb = A * P * A' + Q;     %P(k+1)
    J =  H * Pb * H' + R;    %holder variable R = 1
    K = Pb * H' * inv(J);    %kalman gainzzzz
    P = Pb - K * H * Pb;     %Post measurement P
    xk(:,i+1) = A*xkm1 + K*(zest(:,i)-H*A*xkm1);
    xkm1 = xk(:,i+1);
%    Pcheck(i+1) = P;
 end
%% check KF in lifted state
xkback = wrapToPi(real(ifft(xk)));
xliftback = wrapToPi(real(ifft(Ylift)));
for i = 1:size(A)
    for k = 1:o2
        errstate2(i,k) =  xkback(i,k) - xliftback(i,k);
    end
end


 figure
for i = 1:size(A)
%     if i == 10
%         plot(errstate2(i,:),'+')
%     elseif i == 11
%         plot(errstate2(i,:),'v') 
%     else
        plot(errstate2(i,:))
%     end
    hold on 
end
legend
figure
plot(xliftback(1,:))
hold on
plot(xliftback(2,:))
hold on
plot(xkback(1,:))
hold on
plot(xkback(2,:))
legend('\theta lift','\theta'' lift','\theta filtered','\theta'' filtered')

%% Run DMD again 
%[Phi ,omega ,lambda ,b,Xdmd] = DMD(xk,14,dt);
X1 = xk (:,1:end -1);
X2 = xk (:,2:end);
%% DMD
[U, S, V] = svd(X1, 'econ');
r = min(r, size(U,2));
U_r = U(:, 1:r); % truncate to rank -r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * X2 * V_r / S_r; % low-rank dynamics
[W_r , D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
x1 = X1(:, 1);
% b = Phi\ifft2(x1,order,1);
b = Phi\x1
mm1 = size(X1, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
Xdmdf = Phi * time_dynamics ;

figure
plot (yn(1,1:length(tspan)-1))
hold on 
plot (yn(2,1:length(tspan)-1))
hold on 
xfinal = wrapToPi(real(ifft(Xdmdf)));
plot(xfinal(1,1:length(tspan)-1))
hold on 
plot(xfinal(2,1:length(tspan)-1))
%%
%Bring back lifted coordinates to initial state/

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
%  DataStore.Augx = augx;
 DataStore.Atilda = A;
%  DataStore.Btilda = B2;
 %DataStore.XDMD = Xdmd;
% DataStore.XDMDC = Xdmdc;
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

error = y1check(end,:);
erroravg = [mean(y1check(:,1)) mean(y1check(:,2)) mean(y1check(:,3)) mean(y1check(:,4))];
end

