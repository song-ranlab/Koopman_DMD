clear all; close all; clc
%This is just a function Test Bed

dt = 0.001;
duration = 20;
tspan = 0.0:dt:duration;

% %% Quadrotor DMD Test
% x0 = [ 5;5;5;0;0;0;0;0;0;0;0;0];
% xf = [0;0;0;0;0;0;0;0;0;0;0;0];
%  name = DMDquadrotor(x0,xf,duration
%   

%% DMDc Algorithm test
% 
% X =     [   4 2 1 .5 .25;
%             7 .7 .07 .007 .0007];
%     
% gamma = [ 0 -4 -2 -1 -.5 ];
%  B = [1;0];
% [Mode2,ceval2,deval2,magmode2,Xdmd2,Atild,Xp] = DMDcTest(X,gamma,B,dt);

%% DMDc Algorith test even though nonlinear
% 
% % x0 = [ -3;0;1.01;0];
% % xf = [ 1;0;pi;0];
% x0 = [ -3;0;1.01;0];
% xf = [ 1;0;pi;0];
% 
% count = 1;
% [name,enderr,avgerr] = DMDcartpendEDMD(x0,xf,duration,count);

%% EDMD With Kalman Filter

% x0 = [ -3;0;1.01;0];
% xf = [ 1;0;pi;0];
x0 = [ 1/4;0];
xf = [ 0;0];

count = 1;
[name,enderr,avgerr] = DMDpendEDMD(x0,xf,duration,count);


%%
% g = 9.81;
% m = .028;  %Vehicle Mass
% b = 0.01; %base width
% l=.025;%Arm Length
% I_x = 2.3951 *10^-5; %from paper
% I_y = 2.3951 *10^-5;
% I_z = 1.8580 *10^-5;
% t= 1;
% Stat = [I_x; I_y; I_z; m];
% quadtrig(t,x0,[0;0;0;0],Stat)
