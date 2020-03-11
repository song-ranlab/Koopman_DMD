clear all, close all, clc

ModelName = 'Pendulum_';
ModelName1 = [ModelName, 'Uncontrol_m1_M4_0_0_p25pi_0_v2_'];
ModelName2 = [ModelName, 'Uncontrol_Hamiltonian_m1_M4_0_0_p25pi_0_v2_'];
path2data = ['../Data/',ModelName1]; mkdir(path2data)

path2figs = ['../Figures/PENDULUM/',ModelName1,'/']; mkdir(path2figs)
path2figs2 = ['../Figures/PENDULUM/',ModelName2,'/']; %mkdir(path2figs)
% myVideo = VideoWriter([ModelName1,'_video.mp4','MPEG-4']); %open video file
% myVideo.FrameRate = 5;
% myVideo.Quality = 100;
% ModelName1 = [ModelName, 'Uncontrol_V77_'];
g = 9.81;
m = 1;  %Pendulum Mass
M = 4;  %Cart Mass
L = 1.5;%Arm Length
C = 0;  %Damping - Currently Unused
duration = 4;
dt = 0.05;
% Parameters
tspan = 0.0:dt:duration;
x0 = [0;0;.25*pi;1e-6]; % Initial Conditions

Q  = 1;     %Process Noise Covariance (Go Back To) % 20:0.01 Ratio brings the phase plot to symmetric form
R  = 1;    %Signal Noise Covariance 

%% Plot Settings
ivp=[x0(3); x0(4); g; M; L ; C];
arrow = 'ShowArrow';
interval = [0, duration];
% movie = true;
%%

% H = @(x)(.5*(m+M)*x(2)^2 + m*((.5*L^2*x(4)^2)-L*x(2)*x(4)*sin(x(3))+
% g*L*cos(x(3)))); no
% H = @(x)(.5*(4/3*m*L^2)*x(4)+m*L*g*(cos(x(3))-1)); no
H = @(x)(.5* M * x(2)^2 +...
        .5 * m*((x(2)^2 + ...
        2*x(2)*L*x(4)*sin(x(3))+ ...
        (L*x(4)*sin(x(3)))^2)+...
        (L*x(4)*cos(x(3)))^2) +...
        m*g*L*(cos(x(3))+1));
%     H = @(x)(.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1));
% H = @(x)(
% gradH = @(x)([ 0;... %previous +2
%                 M*x(2) + m*x(2) - m*L*x(4)*cos(x(3));...
%                 m/2*(-2*L*x(2)*x(4)^2*sin(x(3))+2*L^2*x(4)^3*sin(x(3))*cos(x(3)))+m*g*L*x(4)*sin(x(3));...
%                 m/2*(2*L*x(2)*cos(x(3))+2*L^2*x(4)*cos(x(3)))+2*L^2*x(4)*sin(x(3))^2]);gradH = @(x)([ 0;...
% gradH = @(x)([ 0;...  %previous
%                 M*x(2) + m*x(2) - m*L*x(4)*sin(x(3));...
%                 -m*g*L*sin(x(3))-m*L*x(2)*x(4)*cos(x(3));...
%                 m*L^2*x(4)-m*L*x(2)*cos(x(3))]);
gradH = @(x)([  0;...
                (M+m)*x(2)+m*L*x(4)*sin(x(3));
                m*L*x(2)*x(4)*cos(x(3))-m*g*L*sin(x(3));
                m*L*x(2)*sin(x(3))+m*L^2*x(4)]);

[U,V] = meshgrid([-2*pi:0.01:2*pi], [-2*pi:0.01:2*pi]);
% [X,Y] = meshgrid([-4:0.01:4], [-2*pi:0.01:2*pi]); %Legacy
% [X,Y] = meshgrid([-4:0.01:4],[-2*pi:0.01:2*pi] );
[Y,X] = meshgrid([-2*pi:0.01:2*pi],[-4:0.01:4] );
% xset = 0:0.1:.8;
% tset = pi*(0:0.1:2);
%% Unforced
REF = 0;
B = [0 ;0; 0;0];
% f = @(t,x,u)([  x(2);...
%                 (L*x(4)^2-g*sin(x(3)))/cos(x(3));...
%                 x(4);...
%                 ((m+M)*g*sin(x(3))-m*L*x(4)^2*sin(x(3))*cos(x(3)))/(L*(M+m*(1-cos(x(3))^2)))
% ]+B*u);
% f = @(t,x,u)([  x(2);...                   %works
%                 (1/((M+m*(1-cos(x(3))^2))))*(m*sin(x(3))*(L*x(4)^2-g*cos(x(3))));...
%                 x(4);...
%                 (g*sin(x(3))*(m + M))/(L*(m + M - m*cos(x(3))^2))...
%                 ]+B*u);
% f = @(t,x,u)([  x(2);...
%                (1/((M+m*(1-cos(x(3))^2))))*(m*sin(x(3))*(L*x(4)^2-g*cos(x(3))));...
%                 x(4);...
%                 (1/(L*(M+m*(1-cos(x(3))^2)))*((m+M)*g*sin(x(3)) - m*L*x(4)^2*cos(x(3)*sin(x(3)))))]+B*u);
% Sx = @(x)(sin(x(3)));
% Cx = @(x)(cos(x(3)));
d=0;
% D = @(x)(m*L*L*(M+m*(1-Cx^2)));
f = @(t,x,u) ([x(2);
               ((1/(m*L*L*(M+m*(1-cos(x(3))^2))))*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + 0*m*L*L*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))*u);
           x(4);
           (1/(m*L*L*(M+m*(1-cos(x(3))^2))))*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - 0*m*L*cos(x(3))*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))*u]+B*u);
               
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
[t,y0] = ode45(@(t,x)f(t,x,0),tspan,x0,ode_options);
[Hvals0,Jvals0] = evalCostFun_Hamiltonian(H,y0,zeros(1,size(y0,1)),Q,R,REF);


%% Controlled system
%  B = [1 ;0; 0;0];
B = [0; 1;0; 0];
% B = [0 ;0; 1;0];
% B = [0 ;0; 0;1];
% ModelName1 = [ModelName, 'B01_'];
% B = [1; 0];

%f = @(t,x,u)([x(2); -sin(x(1))]+B*u);
% 
% 
% KRONIC: Case 1
REF = 20 + H([0;0;0;0]);
Href = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-REF);
% gain = @(x)(lqr(f(0,x,0),(gradH(x)'*B),Q,R))
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*Href(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,REF);
% 
% % % Store results
DataStore.y1 = y1;
DataStore.u1 = uvals1;
DataStore.H1 = Hvals1;
DataStore.J1 = Jvals1;
DataStore.tspan1 = tspan;
% 
% KRONIC: Case 2
REF = H([0;0;0;0]);
Href = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y2] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals2 = zeros(1,length(y2)); 
for k=1:length(y2), uvals2(1,k) = - gain(y2(k,:))*Href(y2(k,:)); end
[Hvals2,Jvals2] = evalCostFun_Hamiltonian(H,y2,uvals2,Q,R,REF);

% Store results
DataStore.y2 = y2;
DataStore.u2 = uvals2;
DataStore.H2 = Hvals2;
DataStore.J2 = Jvals2;
DataStore.tspan2 = tspan;

% % % % KRONIC: Case 3
REF = 20;
Href = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y3] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals3 = zeros(1,length(y3)); for k=1:length(y3), uvals3(1,k) = - gain(y3(k,:))*Href(y3(k,:)); end
[Hvals3,Jvals3] = evalCostFun_Hamiltonian(H,y3,uvals3,Q,R,REF);

Htest = zeros(size(Hvals1));
for k=1:length(y3), Htest(k) = H([mod(y3(k,1),2*pi), y3(k,2),y3(k,3),y3(k,4)]); end % seems ok, =2

% Store results
DataStore.y3 = y3;
DataStore.u3 = uvals3;
DataStore.H3 = Hvals3;
DataStore.J3 = Jvals3;
DataStore.tspan3 = tspan;

% % KRONIC: Case 4
REF = H([0,0,pi,0]); %REF = -1;
Href = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y4] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals4 = zeros(1,length(y4)); for k=1:length(y4), uvals4(1,k) =-gain(y4(k,:))*Href(y4(k,:)); end
[Hvals4,Jvals4] = evalCostFun_KoopEfun(H,y4,uvals4,Q,R,REF);

% Store results
DataStore.y4 = y4;
DataStore.u4 = uvals4;
DataStore.H4 = Hvals4;
DataStore.J4 = Jvals4;
DataStore.tspan4 = tspan;

%% SAVE RESULTS
save([path2data,'/',[ModelName1,'Data.mat']])