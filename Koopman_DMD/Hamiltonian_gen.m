clear all, close all, clc
path2data = '../Data/';
ModelName1='Pendulum_Uncontrol_m1_M5_0_0_p75pi_0_v3_';
load([path2data,'/',ModelName1,'/',[ModelName1,'Data.mat']])
path2figs = ['../Data/',ModelName1,'/Figures/PENDULUM_Hamiltonian/']; mkdir(path2figs)

myVideo = VideoWriter([ModelName1,'Hamiltonian_video.mp4','MPEG-4']); %open video file
myVideo.FrameRate = 5;
myVideo.Quality = 100;
dxmin = min(y0(:,2));
dxmax = max(y0(:,2));
dtmin = min(y0(:,4));
dtmax = max(y0(:,4));
xset = [dxmin:((dxmax-dxmin)/20):dxmax];
tset = [dtmin:((dtmax-dtmin)/20):dtmax];
figure; hold on
open(myVideo)
for i = 1:length(xset)
delete(findall(gcf,'type','annotation'))
    subplot(2,1,1)
    hold on   
   
    Hfield = (1/2 * M * xset(i)^2 + .5 * m*((xset(i)+L*V.*sin(U)).^2+(L*V.*cos(U)).^2) +m*g*L.*(1-cos(U))); % Legacy Code, Currently Unecessary 
%     Hfield = (.5*(m+M)*xset(i)^2 + m*((.5*L^2*V^2)-L*xset(i)*V*sin(U)+ g*L*cos(U)));
    imagesc(U(1,:),V(:,1),Hfield), shading interp, view(2), colormap(gray(10))
%     plot(y0(:,3),y0(:,4),'-','Color','b','LineWidth',1);
    hold on , grid on 
    plot(DataStore.y1(1:end,3),DataStore.y1(1:end,4),'-','Color','m','LineWidth',1.2); 
%     plot(DataStore.y2(1:end,3),DataStore.y2(1:end,4),'-','Color','g','LineWidth',1.2); 
%     plot(DataStore.y3(1:end,3),DataStore.y3(1:end,4),'-','Color','y','LineWidth',1.2);
%     plot(DataStore.y4(1:end,3),DataStore.y4(1:end,4),'-','Color','c','LineWidth',1.2);
    xlabel('\theta'), ylabel('\theta''')
    hold on

    % Separatrices
%     plot(2*atan(sinh([-pi:0.1:pi])), 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) % center
%     plot(-2*atan(sinh([-pi:0.1:pi])), -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
%     plot(2*atan(sinh([-pi:0.1:pi]))-2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %left
%     plot(-2*atan(sinh([-pi:0.1:pi]))-2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
%     plot(2*atan(sinh([-pi:0.1:pi]))+2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %right
%     plot(-2*atan(sinh([-pi:0.1:pi]))+2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
    axis([-2*pi-0.02 2*pi+0.02 -2*pi-0.02 2*pi+0.02])

    set(gca,'xtick',[-2*pi,-pi,0,pi,2*pi],'xticklabel',{'-2\pi', '-\pi', '0', '\pi', '2\pi'})
    set(gca,'ytick',[-2*pi,-pi,0,pi,2*pi],'yticklabel',{'-2\pi', '-\pi', '0', '\pi', '2\pi'})
    set(gca,'FontSize',14, 'LineWidth',1)
    set(gcf,'Position',[0 0 1080 980])
    set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'nextplot','replacechildren');
    %  grid on ;
    set(gcf,'PaperPositionMode','auto')
    str = sprintf('x'' = %f',xset(i));
    dim = [.2 .625 .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    colorbar%('AxisLocation','in')annotation('textbox',dim,'String',str,'FitBoxToText','on');
    % print('-depsc2', '-loose', [path2figs,ModelName,'PhasePlot','.eps']);
    % print('-dpng', '-loose', [path2figs,ModelName,'PhasePlot','.png']);
    % delay(200);
    legend
    drawnow
    %     F(i) = getframe(gcf); 202
    %     saveas(gcf,sprintf('FIG%d.png',i))
    hold off;
    % end
    %%
    % figure;
    % for i = 1:length(tset)
    subplot(2,1,2)
         hold on
     Hfield2 = (1/2 * M .* X.^2 + .5 .* (m*(X+L*tset(i).*sin(Y)).^2+(L.*tset(i).*cos(Y)).^2) +g*L.*+(1-cos(Y)));
%     Hfield2 =(.5*(m+M)*Y.^2 + m*((.5*L^2.*tset(i).^2)-L*Y.*tset(i).*sin(X)+ g*L*cos(X)));
     imagesc(Y(1,:),X(:,1),Hfield2), shading interp, view(2), colormap(gray(10))
    plot(y0(:,3),y0(:,2),'-','Color','b','LineWidth',1);
    hold on , grid on 
%      plot(DataStore.y1(1:end,3),DataStore.y1(1:end,2),'-','Color','g','LineWidth',1.2); 
%      plot(DataStore.y2(1:end,3),DataStore.y2(1:end,2),'-','Color','m','LineWidth',1.2); 
%      plot(DataStore.y3(1:end,3),DataStore.y3(1:end,2),'-','Color','y','LineWidth',1.2);
%      plot(DataStore.y4(1:end,3),DataStore.y4(1:end,2),'-','Color','c','LineWidth',1.2);
    xlabel('\theta'), ylabel('x''')
    hold on

    %Separatrices
%     plot(2*atan(sinh([-pi:0.1:pi])), 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) % center
%     plot(-2*atan(sinh([-pi:0.1:pi])), -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
%     plot(2*atan(sinh([-pi:0.1:pi]))-2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %left
%     plot(-2*atan(sinh([-pi:0.1:pi]))-2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
%     plot(2*atan(sinh([-pi:0.1:pi]))+2*pi, 2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1) %right
%     plot(-2*atan(sinh([-pi:0.1:pi]))+2*pi, -2*sech([-pi:0.1:pi]),'--','Color',0.7*ones(1,3), 'LineWidth',1)
    axis([-2*pi-0.02 2*pi+0.02 -4-0.03 4+0.03])
    % axis([-2*pi-0.02 2*pi+0.02 -2*pi-0.02 2*pi+0.02])

    set(gca,'ytick',[-4 -2,-1,0,1,2 4],'yticklabel',{'-4','-2', '-1', '0', '1', '2','4'})
    set(gca,'xtick',[-2*pi,-pi,0,pi,2*pi],'xticklabel',{'-2\pi', '-\pi', '0', '\pi', '2\pi'})

    set(gca,'FontSize',14, 'LineWidth',1)
    % set(gcf,'Position',[100 100 350 250])
    set(gcf,'PaperPositionMode','auto')
    colorbar%('AxisLocation','in')
    str = sprintf('\theta'' = %f',tset(i));
    dim = [.2 .151 .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    drawnow
    legend
    %     F(i) = getframe(gcf); 
    %     saveas(gcf,sprintf('FIG%d.png',i+length(xset)))
    hold off
        F(i) = getframe(gcf); 
%     saveas(gcf,sprintf('FIG%d.png',i))
    hmp = sprintf('Hamiltonain_FIG%d.png',i);
    saveas(gcf,[path2figs,hmp]);
end
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(myVideo, frame);
end
close(myVideo)