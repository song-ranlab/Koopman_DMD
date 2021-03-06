clear all, close all, clc

ModelName = 'Quadrotor_';
ModelName1 = [ModelName, 'Uncontrol_z10_u0_v0_w0_p0_q0_r0_'];
ModelName2 = [ModelName, 'Uncontrol_Hamiltonian_z10_u0_v0_w0_p0_q0_r0_'];
path2data = ['../Data/',ModelName1]; mkdir(path2data)

path2figs = ['../Figures/QUADROTOR/',ModelName1,'/']; mkdir(path2figs)
path2figs2 = ['../Figures/QUADROTOR/',ModelName2,'/']; %mkdir(path2figs)

g = 9.81;
m = .28;  %Vehicle Mass
b = 0.01; %base width
l=.025;%Arm Length
I_x = 2.3951 *10^-5; %from paper
I_y = 2.3951 *10^-5;
I_z = 1.8580 *10^-5;

C = 0;  %Damping - Currently Unused
duration = 5;
dt = 0.0001;
% Parameters
tspan = 0.0:dt:duration;
x0 = [0;0;-30;...   %x y z remember negative z is up
    1e-10;1e-10;-5;...  
    0;0;0;...
    pi/4;pi/4;1e-10]; % Initial Conditions x y z u v w psi theta phi p q r

Q  = 1; %for LQR     
R  = 1; %for LQR     
interval = [0, duration];
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);%%

% B = [0;0;0;(1/m);0;(1/m);0;0;0;(1/Ix);(1/Iy);(1/Iz)];
% H = @(x,p,u)(   p(1)*(x(6)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) - x(5)*(cos(x(7))*sin(x(9)) - cos(x(9))*sin(x(7))*sin(x(8))) + x(4)*cos(x(8))*cos(x(9))) +...
%                 p(2)*(x(5)*(cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(8))*sin(x(9))) - x(6)*(cos(x(9))*sin(x(7)) - cos(x(9))*sin(x(7))*sin(x(8))) + x(4)*cos(x(8))*sin(x(9))) +...
%                 p(3)*(x(6)*cos(x(7))*cos(x(8)) - x(4)*sin(x(8)) + x(5)*cos(x(8))*sin(x(7))) +...
%                 -(p(4)*u(1)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))))/m +...
%                 p(5)*u(1)*(cos(x(9))*sin(x(7)) - cos(x(7))*sin(x(8))*sin(x(9)))/m +...
%                 p(6)*(g - (u(1)*cos(x(8))*cos(x(9)))/m) +...
%                 p(7)*(x(10) + x(12)*cos(x(7))*tan(x(8)) + x(11)*sin(x(7))*tan(x(8))) +...
%                 p(8)*(x(11)*cos(x(7)) - x(12)*sin(x(7))) +...
%                 p(9)*((x(12)*cos(x(7)))/cos(x(8)) + (x(11)*sin(x(7)))/cos(x(8)) +...
%                 p(10)*(u(2)/I_x + (x(11)*x(12)*(I_y - I_z))/I_x) +...
%                 p(11)*(u(3)/I_y - (x(10)*x(12)*(I_x - I_z))/I_y) +...
%                 p(12)*(u(4)/I_z + (x(10)*x(11)*(I_x - I_y))/I_z) + u(1)^2/2 + u(2)^2/2 + u(3)^2/2 + u(4)^2/2); 


H = @(x) (m/2 * (x(4)^2+x(5)^2+x(6)^2) - m*g*x(3) + 1/2 *(I_x*x(10)^2+I_y*x(11)^2+I_z*x(12)^2));
gradH = @(x)([0;...
                0;...
                0;...
                x(4);
                x(5);
                x(6);
                0;
                0;
                0;
                x(10);
                x(11);
                x(12)]);
% gradH = @(x)( [ diff(H,x(1));
%                 diff(H,x(2));
%                 diff(H,x(3));
%                 diff(H,x(5));
%                 diff(H,x(6));
%                 diff(H,x(7));
%                 diff(H,x(8));
%                 diff(H,x(9));
%                 diff(H,x(10));
%                 diff(H,x(11));
%                 diff(H,x(12))]);                


%% Unforced
H_ref = 0;
% B = [0 ;0; 0;0];

Fw = 0;   %Wind Force
Tw = 0;   %Wind Torque
Ft = 0;   %Thrust Force
Tx = 0;   %Tau_x
Ty=0;     %Tau_y
Tz=0;     %Tau_z
% Ft = b * (w1^2 + w2^2 + w3^3 + w4^2); %Thrust Force
% Tx = b * l * (w3^2 - w1^2); % X axis torque or roll
% Ty = b * l * (w4^2 - w2^2); % Y axis torque or pitch
% Tz = d * (w2^2 + w4^2 - w1^2 - w3^2); % Z axis torque or yaw
% B = [0;0;0;(1/m);0;(1/m);0;0;0;(1/I_x);(1/I_y);(1/I_z)];
% f = @(t,x,u) ([ x(4);
%                 x(5);
%                 x(6);
%                 (-g*sin(x(8)) + x(12)*x(5) - x(11)*x(6) + 0*(Ft)/m);
%                 (g*cos(x(8))*sin(x(7)) + x(10)*x(8) - x(12)*x(4) + Fw/m);
%                 (g*cos(x(8))*cos(x(7)) - x(10)*x(5) + x(11)*x(4) + Fw/m -0*(Ft)/m);
%                 x(10);
%                 x(11);
%                 x(12);
%                 ((0*(Tx) + Tw + x(11)*x(12) * (Iy-Iz)) / Ix);
%                 ((0*(Ty) + Tw + x(10)*x(12) * (Iz-Ix)) / Iy);
%                 ((0*(Tz) + Tw + x(10)*x(11) * (Ix-Iy)) / Iz) ]+B*u);
%  f = @(t,x,u) ([ x(6) * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) + x(5) * (sin(x(7))*cos(x(9))*sin(x(8)) - sin(x(9))*cos(x(7))) + x(4) * (cos(x(8))*cos(x(9)));
%                  x(5) * (cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(9))*sin(x(8))) + x(6) * (sin(x(7))*cos(x(9))*sin(x(8)) - cos(x(9))*sin(x(7))) + x(4) * (cos(x(8))*sin(x(9)));
%                  x(6) * (cos(x(7))*cos(x(8))) - x(4) * sin(x(8)) + x(5) * (cos(x(8))*sin(x(7))); 
%                  -u(1)/m * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8)));
%                  -u(1)/m * (cos(x(7))*sin(x(9))*sin(x(8)) - cos(x(9))*sin(x(7)));
%                  g - u(1)/m * (cos(x(8))*cos(x(9)));
%                  x(10) + x(12) * (cos(x(7))*tan(x(8))) + x(11) * (sin(x(7))*tan(x(8)));
%                  x(11)*cos(x(7)) - x(12)*sin(x(7));
%                  x(12)*(cos(x(7))/cos(x(8))) + x(11)*(sin(x(7))/cos(x(8)));
%                  ((I_y-I_z)/I_x) * x(11)*x(12) + u(2)/I_x;
%                  ((I_z-I_x)/I_y) * x(10)*x(12) + u(3)/I_y;
%                  ((I_x-I_y)/I_z) * x(10)*x(11) + u(4)/I_z]+B*u);
g1x = @(x)(-1/m * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))));
g2x = @(x)(-1/m * (sin(x(7))*cos(x(9)) - cos(x(7))*sin(x(9))*sin(x(8))));
g3x = @(x)(-1/m * (cos(x(7))*cos(x(8))));
g1 = @(x)([0;0;0;g1x(x);g2x(x);g3x(x);0;0;0;0;0;0]);
g2 = [0;0;0;0;0;0;0;0;0;(1/I_x);0;0];
g3 = [0;0;0;0;0;0;0;0;0;0;(1/I_y);0];
g4 = [0;0;0;0;0;0;0;0;0;0;0;(1/I_z)];

% % B = @(x)(   [0 0 0 0;
%              0 0 0 0;
%              0 0 0 0;
%              g1x(x) 0 0 0;
%              g2x(x) 0 0 0;
%              g3x(x) 0 0 0;
%              0 0 0 0;
%              0 0 0 0; 
%              0 0 0 0;
%              0 (1/Ix) 0 0; 
%              0 0 (1/Iy) 0; 
%              0 0 0 (1/Iz)]);


f = @(t,x,u) ([ x(4);
                x(5);
                x(6);
                0;
                0;
                g;
                (x(10) + x(12) * (cos(x(7))*tan(x(8))) + x(11) * (sin(x(7))*tan(x(8))));
                (x(11)*cos(x(7)) - x(12)*sin(x(7)));
                (x(12)*(cos(x(7))/cos(x(8))) + x(11)*(sin(x(7))/cos(x(8))));
                ((x(11)*x(12) * (I_y-I_z)) / I_x);
                ((x(10)*x(12) * (I_z-I_x)) / I_y);
                ((x(10)*x(11) * (I_x-I_y)) / I_z)] + (g1(x)*u(1) + g2*u(2) + g3*u(3) + g4*u(4)));
            %(g1*u(1) + g2*u(2) + g3*u(3) + g4*u(4))
            

[t,y0] = ode45(@(t,x)f(t,x,[0;0;0;0]),tspan,x0,ode_options);
[Hvals0,Jvals0] = evalCostFun_Hamiltonian(H,y0,zeros(1,size(y0,1)),Q,R,H_ref);

for i = 1:1000
    checker(i,:) = f(t(i),y0(i,:),[0,0,0,0]);
end
%% Controlled system

% % KRONIC: Case 1
H_ref = H([0;0;0;0;0;0;0;0;0;0;0;0]);
% H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) - m*g*L*(1-cos(x(3))))-H_ref);
H_err = @(x)(H(x)-H_ref);
% gain = @(x)(lqr(f(0,x,0),(gradH(x)'*B),Q,R))
% checker = H_err(x0)
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[t,y1] = ode45(@(t,x)f(t,x,-gain(x)*H_err(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*H_err(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,H_ref);
% 
% Store results
DataStore.y1 = y1;
DataStore.u1 = uvals1;
DataStore.H1 = Hvals1;
DataStore.J1 = Jvals1;
DataStore.tspan1 = tspan;
% 
% % KRONIC: Case 2
% H_ref = H([0;0;0;0]);
% H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-H_ref);
% gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
% [~,y2] = ode45(@(t,x)f(t,x,-gain(x)*H_err(x)),tspan,x0,ode_options);
% uvals2 = zeros(1,length(y2)); 
% for k=1:length(y2), uvals2(1,k) = - gain(y2(k,:))*H_err(y2(k,:)); end
% [Hvals2,Jvals2] = evalCostFun_Hamiltonian(H,y2,uvals2,Q,R,H_ref);
% 
% % Store results
% DataStore.y2 = y2;
% DataStore.u2 = uvals2;
% DataStore.H2 = Hvals2;
% DataStore.J2 = Jvals2;
% DataStore.tspan2 = tspan;
% 
% % KRONIC: Case 3
% H_ref = 20;
% H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-H_ref);
% gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
% [~,y3] = ode45(@(t,x)f(t,x,-gain(x)*H_err(x)),tspan,x0,ode_options);
% uvals3 = zeros(1,length(y3)); for k=1:length(y3), uvals3(1,k) = - gain(y3(k,:))*H_err(y3(k,:)); end
% [Hvals3,Jvals3] = evalCostFun_Hamiltonian(H,y3,uvals3,Q,R,H_ref);
% 
% Htest = zeros(size(Hvals1));
% for k=1:length(y3), Htest(k) = H([mod(y3(k,1),2*pi), y3(k,2),y3(k,3),y3(k,4)]); end % seems ok, =2
% 
% % Store results
% DataStore.y3 = y3;
% DataStore.u3 = uvals3;
% DataStore.H3 = Hvals3;
% DataStore.J3 = Jvals3;
% DataStore.tspan3 = tspan;
% 
% % KRONIC: Case 4
% H_ref = H([0,0,pi,0])%REF = -1;
% H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) + m*g*L*(cos(x(3))+1))-H_ref);
% gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
% [t,y4] = ode45(@(t,x)f(t,x,-gain(x)*H_err(x)),tspan,x0,ode_options);
% uvals4 = zeros(1,length(y4)); for k=1:length(y4), uvals4(1,k) =-gain(y4(k,:))*H_err(y4(k,:)); end
% [Hvals4,Jvals4] = evalCostFun_KoopEfun(H,y4,uvals4,Q,R,H_ref);
% 
% % Store results
% DataStore.y4 = y4;
% DataStore.u4 = uvals4;
% DataStore.H4 = Hvals4;
% DataStore.J4 = Jvals4;
% DataStore.tspan4 = tspan;

%% SAVE RESULTS
save([path2data,'/',[ModelName1,'Data.mat']])


%% Kronic Gen for Unforced

y = y0;
usesine = 0;    %not cetain what this is for
polyorder = 4;
nvar = 4;       %nvar is for variation in n rank was 2 in original code

Hy = zeros(length(y),1);
dy = zeros(size(y));
for k=1:length(y)
    dy(k,:) = f(0,y(k,:),[0;0;0;0])';
    Hy(k) = H(y(k,:));
end
figure; hold on, box on
plot(t,y(:,3),'-','Color',[0,0,0.7],'LineWidth',2)
plot(t,y(:,6),'-','Color',[0,0.7,0],'LineWidth',2)
plot(t,y(:,8),'-','Color',[0,0,1,.5],'LineWidth',2)
plot(t,y(:,11),'-','Color',[0,0.5,1],'LineWidth',2)
legend('Height','Vert_Vel','Pitch','Pitch_Rate')
xlabel('t'), ylabel('xi')
set(gca,'xtick',[0:2:10])
set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
title('State Plot')
figure
plot3(y(:,1),y(:,2),-y(:,3))
xlabel('X') ;ylabel('Y') ;zlabel('Z')
title('3-D Plot of X,Y and Z')
saveas(figure(2),[path2figs, 'plot3.fig'])
% print('-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_Trajectory','.eps']);

% %% Determining Eigenvalues and Eigenfunctions for KRONIC
% % Construct libraries  - I want to change this to a foier based one since
% % polynomial order library is not realistic due to its interaction with
% % trigonometric functions as seen in the Hamiltonian
% Theta = buildTheta(y,nvar,polyorder,usesine);
% Gamma = buildGamma(y,dy,nvar,polyorder,usesine);
% 
% % Compute SVD   -I am not sure wh theta is canceled out, wouldnt this
% % difference be required
% [U,S,V] = svd(1*Theta - Gamma,'econ');
% 
% % Least-squares Koopman
% K = pinv(Theta)*Gamma;    %pinv is Mooroe-Penros Puedoiners
% K(abs(K)<1e-12) = 0;
% [T,D] = eig(K);
% D = diag(D);
% [~,IX] = sort(abs(D),'ascend');
% 
% % Compute eigenfunction
% xi0 = V(:,end);             % from SVD
% xi0(abs(xi0)<1e-12) = 0; %if sufficiently small value is 0
% 
% D(IX(1))
% xi1 = T(:,IX(1));%+Tls(:,IX(2));  % from least-squares fit
% xi1(abs(xi1)<1e-12) = 0; 
% 
% % % Print coefficients % lets ignore tis for now, not entirely certain what
% % % this does
% % poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
% % poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);
% 
% % Plot evolution of eigenfunction = Hamiltonian
% if length(Hy)~=length(t)
%     t = 1:length(Hy);
% end
%% DMD

[Mode,ceval,deval,magmode,Xdmd] = DMD (y0(:,1:end-1),y0(:,2:end),12,dt);        %trying the first 2 nodes 



