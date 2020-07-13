clear all, close all, clc

T = 20;
count = 0;
myVideo = VideoWriter([ModelName1,'Hamiltonian_video.mp4','MPEG-4']); %open video file
myVideo.FrameRate = 5;
myVideo.Quality = 100;
open(myVideo)
for i = -pi/2:pi/12:pi/2
    for k =-pi/2:pi/12:pi/2
        count = count+1;
        x0 = [i;k];
        xf = [pi;0];
        fig(count) = DMDPend(x0,xf,T);
        
    end
end

for i=1:length(count)
    % convert the image to a frame
    frame = fig(i) ;    
    writeVideo(myVideo, frame);
end
close(myVideo)