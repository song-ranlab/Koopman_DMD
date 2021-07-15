function [name] = plotend(name,Duration,animate)
% close all
path2data = '../Data';
load([path2data,'/',name,'/',[name,'Data.mat']])
path2figs = ['../Data/',name,'/INVPENDULUM/']; mkdir(path2figs)

% nframes=duration*fps; % Number of Frames 

%% Setup
t = 0:dt:Duration;
t = double(t);
state = y2;
% Position and Angular velocity
x = state(1,:);
dtx = state(2,:);
phi = wrapToPi(state(3,:));
dtphi = state(4,:);
staten = y1;
xn = staten(:,1)';
dtxn = staten(:,2)';
phin = wrapToPi(staten(:,3)');
dtphin = staten(:,4)';
U1 = u1;
U2 = u2';

minw = 1.5*min(dtx) ; maxw = 1.5*max(dtx);                  %xdot
minx = 1.5*min(x) ; maxx = 1.5*max(x);
minv = 1.1*min(dtphi) ;
maxv = 1.1*max(dtphi);    
minwn = 1.5*min(dtxn) ; maxwn = 1.5*max(dtxn);                  %xdot
minxn = 1.5*min(xn) ; maxxn = 1.5*max(xn);
minu = -4 ;
maxu = 4;                  %theta
minvn = 1.1*min(dtphin) ;
maxvn = 1.1*max(dtphin);   %theta_dot
minh = 1.1*min(U1) ; maxh = 1.1*max(U1)+1;
minj = 1.1*min(U2) ; maxj = 1.1*max(U2)+1;

fh = figure(1) ;
% set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','menubar','none') ;
set(fh,'name','Inverted Pendulum','numbertitle','off','color', 'w','position', [0 0 1920 1080]) ;
stop = uicontrol('style','toggle','string','stop','background','w');

%% Plot Original States
subplot(5,2,3) ;
h2 = plot(t,xn,'LineWidth',1,'Color','b') ;
% axis([gmin duration minx maxx]) ;
xlim([t(1) t(end)])
ylim([minxn maxxn])
xlabel('t') ;ylabel('x') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','b');

% Plot for time Vs. Velocity 
subplot(5,2,5);
h3 = plot(t,dtxn,'LineWidth',1,'Color','r') ;
% axis([gmin duration minw maxw]) ;
xlim([t(1) t(end)])
ylim([minwn maxwn])
xlabel('t') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Velocity Plot','Color','r');

% Plot for time Vs. Angle 
subplot(5,2,7) ;
h4 = plot(t,phin,'LineWidth',1,'Color','b') ;
% axis([gmin duration minu maxu]) ;
xlim([t(1) t(end)])
ylim([minu maxu])
xlabel('t') ;ylabel('\theta') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angle Plot','Color','b');

% Plot for time Vs. Angular Rate 
subplot(5,2,9) ;
h5 = plot(t,dtphin,'LineWidth',1,'Color','k') ;
% axis([gmin duration minv maxv]) ;
xlim([t(1) t(end)])
ylim([minvn maxvn])
xlabel('t') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angular Rate Plot','Color','k');

%%Plot Reconstructed States
% Plot for time Vs. States 
subplot(5,2,4) ;
h6 = plot(t,x,'LineWidth',1,'Color','b') ;
% axis([gmin duration minx maxx]) ;
xlim([t(1) t(end)])
ylim([minx maxx])
xlabel('t') ;ylabel('x') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','b');

% Plot for time Vs. Velocity 
subplot(5,2,6);
h7 = plot(t,dtx,'LineWidth',1,'Color','r') ;
% axis([gmin duration minw maxw]) ;
xlim([t(1) t(end)])
ylim([minw maxw])
xlabel('t') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Velocity Plot','Color','r');

% Plot for time Vs. Angle 
subplot(5,2,8) ;
h8 = plot(t,phi,'LineWidth',1,'Color','b') ;
% axis([gmin duration minu maxu]) ;
xlim([t(1) t(end)])
ylim([minu maxu])
xlabel('t') ;ylabel('\theta') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angle Plot','Color','b');

% Plot for time Vs. Angular Rate 
subplot(5,2,10) ;
h9 = plot(t,dtphi,'LineWidth',1,'Color','k') ;
% axis([gmin duration minv maxv]) ;
xlim([t(1) t(end)])
ylim([minv maxv])
xlabel('t') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angular Rate Plot','Color','k');

% Plot for time Vs. Original Control
subplot(5,2,1) ;
h10 = plot(t,U1,'LineWidth',1,'Color','m') ;
% axis([gmin duration minh maxh]) ;
xlim([t(1) t(end)])
ylim([minh maxh])
xlabel('t') ;ylabel('u_x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Original Control','Color','m');

% Plot for time Vs. Uval 
subplot(5,2,2) ;
h11 = plot(t,U2,'LineWidth',1,'Color','m') ;
% legend
% axis([gmin duration gmin gmax]) ;
xlim([t(1) t(end)])
ylim([minj maxj])
xlabel('t') ;ylabel('u') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Learned Control','Color','m');


%% save figure of states
saveas(figure(1),[path2data, 'States_n_control.fig'])
saveas(figure(1),[path2figs, 'States_n_control.png'])
%% Animatio starts
if (animate)
    myVideo = VideoWriter([name,'_video.mp4'],'MPEG-4'); %open video file
    myVideo.FrameRate = 5;
    myVideo.Quality = 100;
    fh2 = figure(2) ;
    hold on
    open(myVideo)
    for i=1:5:length(tspan)-1

        subplot(2,1,1)
        drawcartpend_bw_mine(y1(i,1),y1(i,3),mn,M,L)
        title('Original State')
        subplot(2,1,2)
        drawcartpend_bw_mine(y2(1,i),y2(3,i),mn,M,L)    
        title('Reconstructed State')

        writeVideo(myVideo, getframe(gcf));
    end

    close(myVideo)
end

end

