function [name] = errorplot(enderror,avgerror,name)
%UNTITLED2 Summary of this function goes here
close all


datit = sprintf('Inverted pendulum DMD check');
figure('NumberTitle', 'off','visible','on', 'Name', datit);


subplot(4,1,1)
plot(enderror(1,:))
hold on
plot(avgerror(1,:))
hold off
ylim([-1 100])
title('Control Error - X')
legend ('End Error','Avg Error')
ylabel('X');
xlabel('Test Number');

subplot(4,1,2)
plot(enderror(2,:))
hold on
plot(avgerror(2,:))
hold off
ylim([-1 50])
title('Control Error - X''')
legend ('End Error','Avg Error')
ylabel('X''');
xlabel('Test Number');

subplot(4,1,3)
plot(enderror(3,:))
hold on
plot(avgerror(3,:))
hold off
ylim([-1 20])
title('Control Error - \theta')
legend ('End Error','Avg Error')
ylabel('\theta');
xlabel('Test Number');

subplot(4,1,4)
plot(enderror(4,:))
hold on
plot(avgerror(4,:))
hold off
ylim([-1 2])
title('Control Error - \theta''')
legend ('End Error','Avg Error')
ylabel('\theta''');
xlabel('Test Number');

 
saveas(figure(1),[name, '_plot3.fig'])
end

