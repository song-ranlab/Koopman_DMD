function hi = invpendplot(ModelName1)
close all
path2data = '../Data';
load([path2data,'/',ModelName1,'/',[ModelName1,'Data.mat']])

%%
t = 0.0:dt:duration;
state = y1;
% Position and Angular velocity
x = state(:,1)';
dtx = state(:,2)';
phi = state(:,3)';
dtphi = state(:,4)';

% To set the Range os Phase plane, time vs. depl plots
minu = 1.5*min(phi) ;
maxu = 1.5*max(phi);                  %theta
minv = 1.5*min(dtphi) ;
maxv = 1.5*max(dtphi);              %theta_dot
% minu = 0 ; maxu = 1.1*max(abs(phi));
% minv = 0 ; maxv = 1.1*max(abs(dtphi));
minw = 1.5*min(dtx) ; maxw = 1.5*max(dtx);                  %xdot
minx = 1.5*min(x) ; maxx = 1.5*max(x);
gmin = 0; 
fh = figure ;
% set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','menubar','none') ;
set(fh,'name','Inverted Pendulum','numbertitle','off','color', 'w','position', [0 0 1080 1920]) ;
stop = uicontrol('style','toggle','string','stop','background','w');

%% Plot for Phase plane
subplot(4,2,1) ;
h1 = plot(state(:,3),state(:,4),'LineWidth',1,'Color','m') ;
axis([minu maxu minv maxv]) ;
xlabel('\theta') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on 
title('Phase Plane Plot \theta vs \theta''','Color','m')


subplot(4,2,3) ;
h8 = plot(state(:,3),state(:,2),'LineWidth',1,'Color','r') ;
axis([minu maxu minw maxw]) ;
xlabel('\theta') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x''vs \theta','Color','r')

subplot(4,2,5) ;
h9 = plot(state(:,1),state(:,2),'LineWidth',1,'Color','b') ;
axis([minw maxw minv maxv]) ;
xlabel('x') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x vs x''','Color','b')

subplot(4,2,7) ;
h10 = plot(state(:,2),state(:,4),'LineWidth',1,'Color','b') ;
axis([minw maxw minv maxv]) ;
xlabel('x''') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
%set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x'' vs \theta''','Color','b')

%% Plot for time Vs. States 
subplot(6,2,2) ;
h2 = plot(t,state(:,1),'LineWidth',1,'Color','b') ;
% axis([0 duration minx maxx]) ;
xlabel('t') ;ylabel('x') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','c');

% Plot for time Vs. Velocity 
subplot(6,2,4);
h3 = plot(t,state(:,2),'LineWidth',1,'Color','r') ;
% axis([0 duration minw maxw]) ;
xlabel('t') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Velocity Plot','Color','r');

% Plot for time Vs. Angle 
subplot(6,2,6) ;
h4 = plot(t,state(:,3),'LineWidth',1,'Color','b') ;
% axis([0 duration 0 maxu]) ;
xlabel('t') ;ylabel('\theta') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angle Plot','Color','b');

% Plot for time Vs. Angular Rate 
subplot(6,2,8) ;
h5 = plot(t,state(:,4),'LineWidth',1,'Color','g') ;
% axis([0 duration minv maxv]) ;
xlabel('t') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angular Rate Plot','Color','g');

saveas(figure(1),[path2figs, 'Endstateplot.fig'])

% figure;
% 
% plot3(phi,dtphi,dtx)
% xlabel('\theta') ;ylabel('\theta''') ;zlabel('x''')
% title('3-D Plot of \theta, \theta'' and x''')
% saveas(figure(2),[path2figs, 'plot3.fig'])
end

