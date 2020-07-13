function [] = Xerrorplot(X,Xp,error,figtit)

figure('NumberTitle', 'off', 'Name', figtit,'position', [0 0 1080 1920]);

[x,y] = size(X);
for i = 1:x
    
    subplot(x,1,i)
    plot(error(i,:))
    hold on 
    plot(Xp(i,:))
    hold on
    plot(X(i,:))
    % hold on
    % plot(Xp(1,:),'*')
    hold off
    title(['Control Error - X_',num2str(i)])
    legend ('Error','DMD','Normal')
    ylabel(['X_',num2str(i)]);
    xlabel('t(ms)');
end

sgtitle(figtit)
end

