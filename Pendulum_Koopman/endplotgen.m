clear all, close all, clc
path2data = '../Data/';
ModelName1='Pendulum_Uncontrol_m1_M5_0_0_p75pi_0_v3_';
load([path2data,'/',ModelName1,'/',[ModelName1,'Data.mat']])
path2figs = ['../Data/',ModelName1,'/Figures/PENDULUM/']; mkdir(path2figs)

%%
t = 0.0:dt:duration;

% Position and Angular velocity
x = y4(:,1);
dtx = y4(:,2);
phi = y4(:,3)';
dtphi = y4(:,4)';
L = ivp(5); 
H1 = Hvals4;
J1 = Jvals4;
% To set the Range os Phase plane, time vs. depl plots
minu = 1.5*min(phi) ;
maxu = 1.5*max(phi);                  %theta
minv = 1.5*min(dtphi) ;
maxv = 1.5*max(dtphi);              %theta_dot
% minu = 0 ; maxu = 1.1*max(abs(phi));
% minv = 0 ; maxv = 1.1*max(abs(dtphi));
minw = 1.5*min(dtx) ; maxw = 1.5*max(dtx);                  %xdot
minx = 1.5*min(x) ; maxx = 1.5*max(x);
minh = 1.5*min(H1) ; maxh = 1.5*max(H1);
minj = 1.5*min(J1) ; maxj = 1.5*max(J1);
%J2 = Jvals2;
%minj2 = 1.0*min(J2) ; maxj2 = 1.5*max(J2);
%J3 = Jvals3;
%minj3 = 1.0*min(J3) ; maxj3 = 1.5*max(J3);
%J4 = Jvals4;
%minj4 = 1.0*min(J4) ; maxj4 = 1.5*max(J4);
%gmax = max([maxj, maxj2, maxj3, maxj4]);
gmin = 0;
% 
fh = figure ;
% set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','menubar','none') ;
set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','position', [0 0 1080 1920]) ;
stop = uicontrol('style','toggle','string','stop','background','w');



%% Plot for Phase plane
subplot(4,2,1) ;
h1 = plot(y4(:,3),y4(:,4),'LineWidth',1,'Color','m') ;
axis([minu maxu minv maxv]) ;
xlabel('\theta') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on 
title('Phase Plane Plot \theta vs \theta''','Color','m')


subplot(4,2,3) ;
h8 = plot(y4(:,3),y4(:,2),'LineWidth',1,'Color','r') ;
axis([minu maxu minw maxw]) ;
xlabel('\theta') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x''vs \theta','Color','r')

subplot(4,2,5) ;
h9 = plot(y4(:,1),y4(:,2),'LineWidth',1,'Color','b') ;
axis([minw maxw minv maxv]) ;
xlabel('x') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x vs x''','Color','b')

subplot(4,2,7) ;
h10 = plot(y4(:,2),y4(:,4),'LineWidth',1,'Color','b') ;
axis([minw maxw minv maxv]) ;
xlabel('x''') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x'' vs \theta''','Color','b')

%% Plot for time Vs. States 
subplot(6,2,2) ;
h2 = plot(t,y4(:,1),'LineWidth',1,'Color','b') ;
axis([0 duration minx maxx]) ;
xlabel('t') ;ylabel('x') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','c');

% Plot for time Vs. Velocity 
subplot(6,2,4);
h3 = plot(t,y4(:,2),'LineWidth',1,'Color','r') ;
axis([0 duration minw maxw]) ;
xlabel('t') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Velocity Plot','Color','r');

% Plot for time Vs. Angle 
subplot(6,2,6) ;
h4 = plot(t,y4(:,3),'LineWidth',1,'Color','b') ;
axis([0 duration minu maxu]) ;
xlabel('t') ;ylabel('\theta') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angle Plot','Color','b');

% Plot for time Vs. Angular Rate 
subplot(6,2,8) ;
h5 = plot(t,y4(:,4),'LineWidth',1,'Color','g') ;
axis([0 duration minv maxv]) ;
xlabel('t') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angular Rate Plot','Color','g');

% Plot for time Vs. Hval
subplot(6,2,10) ;
h6 = plot(t,H1,'LineWidth',1,'Color','m') ;
% axis([0 duration minh maxh]) ;
xlabel('t') ;ylabel('u_x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Hval','Color','m');

% Plot for time Vs. Uval 
subplot(6,2,12) ;
h7 = plot(t,J1,'LineWidth',1,'Color','k') ;
% hold on 
% h10 = plot(t,J2,'LineWidth',1,'Color','g') ;
% hold on
% h11 = plot(t,J3,'LineWidth',1,'Color','r') ;
% hold on
% h12 = plot(t,J4,'LineWidth',1,'Color','b') ;
% hold off
% legend
axis([0 duration gmin 2000]) ;
xlabel('t') ;ylabel('u_\phi''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs.JVal','Color','k');

saveas(figure(1),[path2figs, 'Endstateplot.fig'])

figure;

plot3(phi,dtphi,dtx)
xlabel('\theta') ;ylabel('\theta''') ;zlabel('x''')
title('3-D Plot of \theta, \theta'' and x''')
saveas(figure(2),[path2figs, 'plot3.fig'])