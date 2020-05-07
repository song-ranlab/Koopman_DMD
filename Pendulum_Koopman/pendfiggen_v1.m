7clear all, close all, clc
path2data = '../Data/';
ModelName1='Pendulum_Uncontrol_m1_M5_0_0_p75pi_0_v2_';
load([path2data,'/',ModelName1,'/',[ModelName1,'Data.mat']])
path2figs = ['../Data/',ModelName1,'/Figures/PENDULUM/']; mkdir(path2figs)

myVideo = VideoWriter([ModelName1,'_video.mp4'],'MPEG-4'); %open video file
myVideo.FrameRate = 5;
myVideo.Quality = 100;

%%
fps = 10;
nframes=duration*fps; % Number of Frames
t = 0.0:dt:duration;

% Position and Angular velocity
x = y0(:,1);
dtx = y0(:,2);
phi = y0(:,3)';
dtphi = y0(:,4)';
L = ivp(5); 
H1 = Hvals0;
J1 = Jvals1;
% To set the Range os Phase plane, time vs. depl plots
minu = 1.1*min(phi) ;
maxu = 1.1*max(phi);                  %theta
minv = 1.1*min(dtphi) ;
maxv = 1.1*max(dtphi);              %theta_dot
% minu = 0 ; maxu = 1.1*max(abs(phi));
% minv = 0 ; maxv = 1.1*max(abs(dtphi));
minw = 1.5*min(dtx) ; maxw = 1.1*max(dtx);                  %xdot
minx = 1.5*min(x) ; maxx = 1.1*max(x);
minh = 1.1*min(H1) ; maxh = 1.1*max(H1);
minj = 1.1*min(J1) ; maxj = 1.1*max(J1);
J2 = Jvals2;
minj2 = 1.1*min(J2) ; maxj2 = 1.1*max(J2);
J3 = Jvals3;
minj3 = 1.1*min(J3) ; maxj3 = 1.1*max(J3);
J4 = Jvals4;
minj4 = 1.1*min(J4) ; maxj4 = 1.1*max(J4);
gmax = max([maxj, maxj2, maxj3, maxj4]);
gmin = min([minj, minj2,minj3,minj4]);
% 
fh = figure ;
% set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','menubar','none') ;
set(fh,'name','The Simple Pendulum','numbertitle','off','color', 'w','position', [0 0 1080 1980]) ;
stop = uicontrol('style','toggle','string','stop','background','w');


% %% Plot for Pendulum
% subplot(4,2,1);
% h = plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',1.5,'Color','b');
% title('Pendulum Animation','Color','b');
% range = 1.1*L;
% axis([-range range -range range]);
% axis square;
% set(gca,'XTickLabelMode', 'manual', 'XTickLabel', [],'YTickLabelMode', .....
%     'manual', 'YTickLabel', [],'nextplot','replacechildren');



%% Plot for Phase plane
subplot(4,2,3) ;
h1 = plot(x(3),x(4),'LineWidth',1,'Color','m') ;
axis([minu maxu minv maxv]) ;
xlabel('\theta') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'nextplot','replacechildren');
grid on 
title('Phase Plane Plot \theta vs \theta''','Color','m')


subplot(4,2,5) ;
h8 = plot(x0(3),x0(2),'LineWidth',1,'Color','r') ;
axis([minu maxu minw maxw]) ;
xlabel('\theta') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x''vs \theta','Color','r')


subplot(4,2,7) ;
h9 = plot(x0(1),x0(2),'LineWidth',1,'Color','b') ;
axis([minw maxw minv maxv]) ;
xlabel('x''') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
set(gca,'nextplot','replacechildren');
grid on ;
title('Phase Plane Plot x'' vs \theta''','Color','b')

%% Plot for time Vs. States 
subplot(6,2,2) ;
h2 = plot(t(1),x(1),'LineWidth',1,'Color','b') ;
axis([0 duration minx maxx]) ;
xlabel('t') ;ylabel('x') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Displacement Plot','Color','c');

% Plot for time Vs. Velocity 
subplot(6,2,4);
h3 = plot(t(1),dtx(1),'LineWidth',1,'Color','r') ;
axis([0 duration minw maxw]) ;
xlabel('t') ;ylabel('x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Velocity Plot','Color','r');

% Plot for time Vs. Angle 
subplot(6,2,6) ;
h4 = plot(t(1),phi(1),'LineWidth',1,'Color','b') ;
axis([0 duration minu maxu]) ;
xlabel('t') ;ylabel('\theta') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angle Plot','Color','b');

% Plot for time Vs. Angular Rate 
subplot(6,2,8) ;
h5 = plot(t(1),dtphi(1),'LineWidth',1,'Color','g') ;
axis([0 duration minv maxv]) ;
xlabel('t') ;ylabel('\theta''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Angular Rate Plot','Color','g');

% Plot for time Vs. Hval
subplot(6,2,10) ;
h6 = plot(t(1),H1,'LineWidth',1,'Color','m') ;
axis([0 duration minh maxh]) ;
xlabel('t') ;ylabel('u_x''') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Hval','Color','m');

% Plot for time Vs. Uval 
subplot(6,2,12) ;
h7 = plot(t(1),J1,'LineWidth',1,'Color','k') ;
% hold on 
h10 = plot(t(1),J2,'LineWidth',1,'Color','g') ;
% hold on
h11 = plot(t(1),J3,'LineWidth',1,'Color','r') ;
% hold on
h12 = plot(t(1),J4,'LineWidth',1,'Color','b') ;
% hold off
% legend
axis([0 duration gmin gmax]) ;
xlabel('t') ;ylabel('Gain') ;
set(get(gca,'YLabel'),'Rotation',0.0)
grid on ;
set(gca,'nextplot','replacechildren');
title('Time Vs. Gain','Color','k');

%% Animatio starts
 
open(myVideo)

for i=1:length(phi)-1
 % Animation Plot
%     if (ishandle(h)==1)
%         Xcoord=[0,L*sin(-phi(i))];
%         Ycoord=[0,L*cos(-phi(i))];
%         set(h,'XData',Xcoord,'YData',Ycoord);
%         if get(stop,'value')==0
%             drawnow;
%         elseif get(stop,'value')==1
%             break;
%         end
%     end
    subplot(4,2,1)
    drawcartpend_bw_mine(y0(i,1),y0(i,3),m,M,L)
    % Phase Plane Plot
    if (ishandle(h1)==1)
        PP(i,:) = [phi(i) dtphi(i)];
        set(h1,'XData',PP(:,1),'YData',PP(:,2));
        drawnow;    
        vasu = length(PP(:,1)) ;
        if arrow == 'ShowArrow'
        if i>1
             subplot(423)
             arrowh(phi(i-1:vasu),dtphi(i-1:vasu),'m') ;
             hold on ;
        end         
        end

        XX(i,:) = [phi(i) dtx(i)];
        set(h8,'XData',XX(:,1),'YData',XX(:,2));
        drawnow;    
        vasu = length(XX(:,1)) ;
        if arrow == 'ShowArrow'
        if i>1
             subplot(425)
             arrowh(phi(i-1:vasu),dtx(i-1:vasu),'r') ;
             hold on ;
        end         
        end
        XT(i,:) = [dtx(i) dtphi(i)];
        set(h9,'XData',XT(:,1),'YData',XT(:,2));
        drawnow;    
        vasu = length(XT(:,1)) ;
        if arrow == 'ShowArrow'
            if i>1
                 subplot(427)
                 arrowh(dtx(i-1:vasu),dtphi(i-1:vasu),'r') ;
                 hold on ;
            end         
        end
    end
    % Time Vs. displacement Plot 
        if (ishandle(h2)==1)
            DEPL(i,:) = [t(i) x(i)] ;
            set(h2,'Xdata',DEPL(:,1),'YData',DEPL(:,2)) ;
            drawnow ;    
        end    
        if (ishandle(h3)==1)
            DEPL1(i,:) = [t(i) dtx(i)] ;
            set(h3,'Xdata',DEPL1(:,1),'YData',DEPL1(:,2)) ;
            drawnow ;    
        end    
        if (ishandle(h4)==1)
            DEPL2(i,:) = [t(i) phi(i)] ;
            set(h4,'Xdata',DEPL2(:,1),'YData',DEPL2(:,2)) ;
            drawnow ;    
        end 
        if (ishandle(h5)==1)
            DEPL3(i,:) = [t(i) dtphi(i)] ;
            set(h5,'Xdata',DEPL3(:,1),'YData',DEPL3(:,2)) ;
            drawnow ;    
        end    
        if (ishandle(h6)==1)
            DEPL4(i,:) = [t(i) H1(i)] ;
            set(h6,'Xdata',DEPL4(:,1),'YData',DEPL4(:,2)) ;
            drawnow ;    
        end    
        if (ishandle(h7)==1)
            DEPL5(i,:) = [t(i) J1(i)] ;
            set(h7,'Xdata',DEPL5(:,1),'YData',DEPL5(:,2)) ;
            drawnow ; 
%             hold on
        end    
%          if (ishandle(h10)==1)
%             DEPL6(i,:) = [t(i) J2(i)] ;
%             set(h12,'Xdata',DEPL6(:,1),'YData',DEPL6(:,2)) ;
%             drawnow ;
% %             hold on
%          end   
%          if (ishandle(h11)==1)
%             DEPL7(i,:) = [t(i) J3(i)] ;
%             set(h11,'Xdata',DEPL7(:,1),'YData',DEPL7(:,2)) ;
%             drawnow ;    
% %             hold on 
%          end   
%          if (ishandle(h12)==1)
%             DEPL8(i,:) = [t(i) J4(i)] ;
%             set(h10,'Xdata',DEPL8(:,1),'YData',DEPL8(:,2)) ;
%             drawnow ;    
% %             hold off
%         end   
        
    
    F(i) = getframe(gcf); 
    %saveas(gcf,sprintf('FIG%d.png',i))
    hmp = sprintf('FIG%d.png',i);
    saveas(gcf,[path2figs,hmp]);
end
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(myVideo, frame);
end
close(myVideo)

figure;

plot3(phi,dtphi,dtx)
xlabel('\theta') ;ylabel('\theta''') ;zlabel('x''')
title('3-D Plot of \theta, \theta'' and x''')
saveas(figure(2),[path2figs, 'plot3.fig'])