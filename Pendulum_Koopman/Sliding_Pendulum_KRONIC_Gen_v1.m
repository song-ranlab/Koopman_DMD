clear all, close all, clc

ModelName = 'Pendulum_';
ModelName1 = [ModelName, 'Uncontrol_m1_M5_0_0_p75pi_0_v3_'];
ModelName2 = [ModelName, 'Uncontrol_Hamiltonian_m1_M5_0_0_p75pi_0_v3_'];
path2data = ['../Data/',ModelName1]; mkdir(path2data)

path2figs = ['../Figures/PENDULUM/',ModelName1,'/']; mkdir(path2figs)
path2figs2 = ['../Figures/PENDULUM/',ModelName2,'/']; %mkdir(path2figs)
% myVideo = VideoWriter([ModelName1,'_video.mp4','MPEG-4']); %open video file
% myVideo.FrameRate = 5;
% myVideo.Quality = 100;
% ModelName1 = [ModelName, 'Uncontrol_V77_'];
g = -9.81;
m = 1;  %Pendulum Mass
M = 5;  %Cart Mass
L = 2;%Arm Length
C = 0;  %Damping - Currently Unused
duration = 10;
dt = 0.001;
% Parameters
tspan = 0.0:dt:duration;
x0 = [-3; 0; pi+.1; 0]; % Initial Conditions
% 
% Q  = 1;     
% R  = 1;  
Q = [1 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100];
R = .0001;  

%% Plot Settings
arrow = 'ShowArrow';
interval = [0, duration];
% movie = true;
%%
% H = @(x)(.5* M * x(2)^2 +...
%     .5 * m*((x(2)^2 + ...
%     2*x(2)*L*x(4)*cos(x(3))+ ...
%     (L*x(4)*sin(x(3)))^2)+...
%     (L*x(4)*cos(x(3)))^2) -...
%     m*g*L*(1-cos(x(3))));
% 
% gradH = @(x)([  0;...
%                 (M+m)*x(2)+m*L*x(4)*sin(x(3));
%                 m*L*x(2)*x(4)*cos(x(3))-m*g*L*sin(x(3));
%                 m*L*x(2)*sin(x(3))+m*L^2*x(4)]);

% [U,V] = meshgrid([-2*pi:0.01:2*pi], [-2*pi:0.01:2*pi]);

% [Y,X] = meshgrid([-2*pi:0.01:2*pi],[-4:0.01:4] );
% xset = 0:0.1:.8;
% tset = pi*(0:0.1:2);
%% Unforced
% H_ref = 0;
% B = [0 ;0; 0;0];
s = 1;
d=1;
A = [0 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0];
B = [0; 1/M; 0; s*1/(M*L)];
f = @(t,x,u)([  x(2);...                   %works
                (1/((M+m*(1-cos(x(3))^2))))*(m*sin(x(3))*(L*x(4)^2-g*cos(x(3))));...
                x(4);...
                (g*sin(x(3))*(m + M))/(L*(m + M - m*cos(x(3))^2))...
                ]+B*u);
% Sx = @(x)(sin(x(3)));
% Cx = @(x)(cos(x(3)));
% d=0;
% D = @(x)(m*L*L*(M+m*(1-Cx^2)));
% f = @(t,x,u) ([x(2);
%                ((1/(m*L*L*(M+m*(1-cos(x(3))^2))))*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + 0*m*L*L*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))*u);
%            x(4);
%            (1/(m*L*L*(M+m*(1-cos(x(3))^2))))*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - 0*m*L*cos(x(3))*(1/(m*L*L*(M+m*(1-cos(x(3))^2))))*u]+B*u);
%                
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
% [t,y0] = ode45(@(t,x)f(t,x,0),tspan,x0,ode_options);
% [Hvals0,Jvals0] = evalCostFun_Hamiltonian(H,y0,zeros(1,size(y0,1)),Q,R,H_ref);


%% Controlled system
%  B = [1 ;0; 0;0];
% B = [0; 1;0; 0];
% B = [0 ;0; 1;0];
% B = [0 ;0; 0;1];
% ModelName1 = [ModelName, 'B01_'];
% B = [1; 0];

%f = @(t,x,u)([x(2); -sin(x(1))]+B*u);
% 
% 
% KRONIC: Case 1
% H_ref = 40 + H([0;0;0;0]);
% H_err = @(x)((.5* M * x(2)^2 + .5 * m*((x(2)^2 + 2*x(2)*L*x(4)*sin(x(3))+ (L*x(4)*sin(x(3)))^2)+ (L*x(4)*cos(x(3)))^2) - m*g*L*(1-cos(x(3))))-H_ref);
%  H_err = @(x)(H(x)-H_ref);
% H_err = @(x)y-[4;0;0;0];
% gain = @(x)(lqr(f(0,x,0),B,Q,R));
% checker = H_err(x0)
gain = lqr(A,B,Q,R);
[~,y1] = ode45(@(t,y1)cartpend(y1,m,M,L,g,d,-gain*(y1-[1; 0; pi; 0])),tspan,x0);
% [~,y1] = ode45(@(t,x)f(t,x,-gain*H_err(x)),tspan,x0,ode_options);
% uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*H_err(y1(k,:)); end
% [Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,H_ref);
% 
% % % Store results
DataStore.y1 = y1;
% DataStore.u1 = uvals1;
% DataStore.H1 = Hvals1;
% DataStore.J1 = Jvals1;
DataStore.tspan1 = tspan;
% % 
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
% % % % % KRONIC: Case 3
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
% % % KRONIC: Case 4
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


% %% Kronic Gen for Unforced
% 
% y = y1;
% usesine = 0;    %not cetain what this is for
% polyorder = 4;
% nvar = 4;       %nvar is for variation in n rank was 2 in original code
% 
% Hy = zeros(length(y),1);
% dy = zeros(size(y));
% for k=1:length(y)
%     dy(k,:) = f(0,y(k,:),0)';
%     Hy(k) = H(y(k,:));
% end
% figure; hold on, box on
% plot(t,y(:,1),'-','Color',[0,0,0.7],'LineWidth',2)
% plot(t,y(:,2),'-','Color',[0,0.7,0],'LineWidth',2)
% plot(t,y(:,3),'-','Color',[0,0,1,.5],'LineWidth',2)
% plot(t,y(:,4),'-','Color',[0,0.5,1],'LineWidth',2)
% legend('x1','x2','x3','x4')
% xlabel('t'), ylabel('xi')
% set(gca,'xtick',[0:2:duration])
% set(gca,'FontSize',16)
% % set(gcf,'Position',[100 100 225 200])
% set(gcf,'PaperPositionMode','auto')
% title('State Plot')
% % print('-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_Trajectory','.eps']);
% 
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
% 
% %% Hamiltonian plot
% figure; hold on, box on
% ph(1) = plot(t,Hy./norm(Hy),'-k', 'LineWidth',18,'Color',[0.7,0.7,0.7]);
% ph(2) = plot(t,(Theta)*(xi0)./norm((Theta)*(xi0)),'-b', 'LineWidth',8,'Color',[0,0,1]);
% ph(3) = plot(t,-(Theta)*(xi1)./norm((Theta)*(xi1)),'--', 'LineWidth',8,'Color',[0,0.5,0]);
% xlabel('t'), ylabel('E')
% ylim([0.01-1e-4 0.2]), xlim([min(t) max(t)])
% set(gca,'xtick',[0:2:duration])
% legend(ph,'True place', 'SVD place', 'LS place')
% title('Hamiltonian Plot')
% set(gca,'FontSize',16)
% % set(gcf,'Position',[100 100 1000 1000])
% set(gcf,'PaperPositionMode','auto')
% % print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);
% 
% % Plot error 
% dstep = 5;
% clear ph
% figure; hold on, box on
% tmp = Gamma*xi0;
% ph(1) = plot(t(1:dstep:end),tmp(1:dstep:end),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
% tmp = -Gamma*xi1;
% ph(2) = plot(t(1:dstep:end),tmp(1:dstep:end),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
% xlabel('t'), ylabel('Gamma xi')
% xlim([min(t) max(t)])
% ylim([min(Gamma*xi1)+0.2*min(Gamma*xi1) max(Gamma*xi1)+0.5*max(Gamma*xi1)])
% legend(ph,'SVD p', 'LS  p')
% set(gca,'xtick',[0:2:10])
% set(gca,'FontSize',16)
% % set(gcf,'Position',[100 100 225 200])
% set(gcf,'PaperPositionMode','auto')
% title('Error Plot')
% % print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);
% 
% err0 = (Hy./norm(Hy)-(Theta)*(xi0)./norm((Theta)*(xi0)))./(Hy./norm(Hy));
% err1 = (Hy./norm(Hy)-(Theta)*(xi1)./norm((Theta)*(xi1)))./(Hy./norm(Hy));
% figure, plot(err0,'-r'),hold on, plot(err1,'--b')
% axis tight
%%
% % return
% 
% %% Ensemble data
% close all
% ModelName_tmp = [ModelName, '_ensemble_'];
% 
% % Collect data
% clear y dy Hy
% [Xgrid,Vgrid] = meshgrid([-2.45:.05:2.5],[-2.45:.05:2.5]);
% Nx = size(Xgrid,1), Ny = size(Vgrid,2);
% y(:,1) = reshape(Xgrid,[Nx*Ny 1]);
% y(:,2) = reshape(Vgrid,[Nx*Ny 1]);
% y(:,3) = reshape(Vgrid,[Nx*Ny 1]);
% y(:,4) = reshape(Vgrid,[Nx*Ny 1]);
% Hy = zeros(length(y),1);
% Hdy = zeros(length(y),1);
% dy = zeros(size(y));
% for k=1:length(y)
%     dy(k,:) = f(0,y(k,:),0)';
%     Hy(k) = H(y(k,:));0
%     Hdy(k) = H(dy(k,:));
% end
% 
% cmap = gray(2*11);
% cmap = cmap(1:2:end,:);
% %xH = [-2.4:0.01:2.4];
% %[X,V] = meshgrid(xH,xH);
% [X,U,V] = meshgrid([-2*pi:0.01:2*pi],[-2*pi:0.01:2*pi], [-2*pi:0.01:2*pi]);
% Hfield = (1/2 * M * X.^2 + .5 * m*((X+L*V.*sin(U)).^2+(L*V.*cos(U)).^2) +m*g*L.*(1-cos(U))) ;
% 
% figure; box on, hold on
% sh = surf(xH,xH,zeros(size(Hfield)));shading interp, view(2)
% sh.FaceColor = [0 0 0];
% contourf(xH,xH,log(Hfield+(0.2500001)),[log([-0.2,-0.1,0,0.25:0.5:4]+(0.2500001))],'LineColor','none'), colormap(cmap)
% hold on
% sh = scatter(Xgrid(:),Vgrid(:),'.b');
% sh.SizeData = 5;
% x1 = -sqrt(2):0.01:sqrt(2);
% plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--y')
% plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--y')
% xlabel('x1'), ylabel('x2')
% drawnow
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 220 200])
% set(gcf,'PaperPositionMode','auto')
% % print('-painters','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);
% % print('-opengl','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);
% 
% % Construct libraries
% Theta = buildTheta(y,nvar,polyorder,usesine);
% Gamma = buildGamma(y,dy,nvar,polyorder,usesine);
% 
% % Compute SVD
% [U,S,V] = svd(0*Theta - Gamma,'econ');
% 
% % Least-squares Koopman
% K = pinv(Theta)*Gamma;
% K(abs(K)<1e-12) = 0;
% [T,D] = eig(K);
% D = diag(D);
% [~,IX] = sort(abs(D),'ascend');
% 
% % Compute eigenfunction
% xi0 = V(:,end);             % from SVD
% xi0(abs(xi0)<1e-12) = 0;
% 
% D(IX(1))
% xi1 = T(:,IX(1));  % from least-squares fit
% xi1(abs(xi1)<1e-12) = 0; 
% 
% % Print coefficients
% poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
% poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);
% 
% % Plot evolution of eigenfunction = Hamiltonian
% if length(Hy)~=length(t)
%     t = 1:length(Hy);
% end
% 
% %% Show results for ensemble of points
% close all
% figure; hold on, box on
% ph(1) = plot(t,Hy./norm(Hy),'-k', 'LineWidth',5,'Color',[0.7,0.7,0.7]);
% ph(2) = plot(t,(Theta)*(xi0)./norm((Theta)*(xi0)),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
% ph(3) = plot(t,-(Theta)*(xi1)./norm((Theta)*(xi1)),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
% xlabel('sampling point'), ylabel('E')
% ylim([-0.005 0.05]),
% xlim([-0 10])
% legend(ph,'True place', 'SVD place', 'LS place','location','north')
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 240 200])
% set(gcf,'PaperPositionMode','auto')
% % print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);
% 
% % Plot error
% dstep = 5;
% clear ph
% figure; hold on, box on
% tmp = Gamma*xi1;
% ph(2) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
% tmp = Gamma*xi0;
% ph(1) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
% xlabel('sampling point'), ylabel('Gamma xi')
% axis tight
% xlim([0,Nx*Ny+10])
% set(gca, 'xtick',[5000 10000],'xticklabel',{'5k', '10k'})
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 225 200])
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);
% 
% return
% %% Ensemble data
% close all
% ModelName_tmp = [ModelName, '_ensemble_'];
% 
% % Collect data
% clear y dy Hy
% [Xgrid,Vgrid] = meshgrid([-2.45:.05:2.5],[-2.45:.05:2.5]);
% Nx = size(Xgrid,1), Ny = size(Vgrid,2);
% y(:,1) = reshape(Xgrid,[Nx*Ny 1]);
% y(:,2) = reshape(Vgrid,[Nx*Ny 1]);
% Hy = zeros(length(y),1);
% Hdy = zeros(length(y),1);
% dy = zeros(size(y));
% for k=1:length(y)
%     dy(k,:) = f(0,y(k,:))';
%     Hy(k) = H(y(k,:));
%     Hdy(k) = H(dy(k,:));
% end
% 
% cmap = gray(2*11);
% cmap = cmap(1:2:end,:);
% xH = [-2.4:0.01:2.4];
% [X,V] = meshgrid(xH,xH);
% Hfield = (1/2)*V.^2-(X.^2)/2 + (1/4)*X.^4; 
% 
% figure; box on, hold on
% sh = surf(xH,xH,zeros(size(Hfield)));shading interp, view(2)
% sh.FaceColor = [0 0 0];
% contourf(xH,xH,log(Hfield+(0.2500001)),[log([-0.2,-0.1,0,0.25:0.5:4]+(0.2500001))],'LineColor','none'), colormap(cmap)
% hold on
% sh = scatter(Xgrid(:),Vgrid(:),'.b');
% sh.SizeData = 5;
% x1 = -sqrt(2):0.01:sqrt(2);
% plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--y')
% plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--y')
% xlabel('x1'), ylabel('x2')
% drawnow
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 220 200])
% set(gcf,'PaperPositionMode','auto')
% % print('-painters','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);
% print('-opengl','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);
% 
% % Construct libraries
% Theta = buildTheta(y,nvar,polyorder,usesine);
% Gamma = buildGamma(y,dy,nvar,polyorder,usesine);
% 
% % Compute SVD
% [U,S,V] = svd(0*Theta - Gamma,'econ');
% 
% % Least-squares Koopman
% K = pinv(Theta)*Gamma;
% K(abs(K)<1e-12) = 0;
% [T,D] = eig(K);
% D = diag(D);
% [~,IX] = sort(abs(D),'ascend');
% 
% % Compute eigenfunction
% xi0 = V(:,end);             % from SVD
% xi0(abs(xi0)<1e-12) = 0;
% 
% D(IX(1))
% xi1 = T(:,IX(1));  % from least-squares fit
% xi1(abs(xi1)<1e-12) = 0; 
% 
% % Print coefficients
% poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
% poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);
% 
% % Plot evolution of eigenfunction = Hamiltonian
% if length(Hy)~=length(t)
%     t = 1:length(Hy);
% end
% 
% %% Show results for ensemble of points
% close all
% figure; hold on, box on
% ph(1) = plot(t,Hy./norm(Hy),'-k', 'LineWidth',5,'Color',[0.7,0.7,0.7]);
% ph(2) = plot(t,(Theta)*(xi0)./norm((Theta)*(xi0)),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
% ph(3) = plot(t,-(Theta)*(xi1)./norm((Theta)*(xi1)),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
% xlabel('sampling point'), ylabel('E')
% ylim([-0.005 0.05]),
% xlim([-0 10])
% legend(ph,'True place', 'SVD place', 'LS place','location','north')
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 240 200])
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);
% 
% % Plot error
% dstep = 5;
% clear ph
% figure; hold on, box on
% tmp = Gamma*xi1;
% ph(2) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
% tmp = Gamma*xi0;
% ph(1) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
% xlabel('sampling point'), ylabel('Gamma xi')
% axis tight
% xlim([0,Nx*Ny+10])
% set(gca, 'xtick',[5000 10000],'xticklabel',{'5k', '10k'})
% set(gca,'FontSize',16)
% set(gcf,'Position',[100 100 225 200])
% set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);
