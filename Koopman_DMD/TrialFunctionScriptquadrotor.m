clear all, close all, clc
ModelName = 'Pendulum_Inverted_';
c = date;
ModelName1 = [ModelName, '_', c, '_'];
T = 30;
count = 0;
myVideo = VideoWriter([ModelName1,'Range_video.mp4','MPEG-4']); %open video file
myVideo.FrameRate = 5;
myVideo.Quality = 100;
open(myVideo)
for i = 1/2:.1:3/2
    for k =-1/4:.1:1/4
        count = count+1
        x0 = [0;0;i;k];
        xf = [5;0;pi;0];
        name = DMDcartpend(x0,xf,T);
        talk = invpendplot(name);
        F(count) = getframe(gcf);
    end
end

for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(myVideo, frame);
end
close(myVideo)