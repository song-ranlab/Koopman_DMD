clear all; close all; clc
%This is just a function Test Bed

dt = 0.001;
duration = 20;
tspan = 0.0:dt:duration;

%% Quadrotor DMD Test
% x0 = [ 5;5;5;0;0;0;0;0;0;0;0;0];
% xf = [0;0;0;0;0;0;0;0;0;0;0;0];
% name = DMDquadrotor(x0,xf,30)

%% DMDc Algorithm test

X =     [   4 2 1 .5 .25;
            7 .7 .07 .007 .0007];
    
gamma = [   -4 -2 -1 -.5];
 B = [0;1];
[Mode2,ceval2,deval2,magmode2,Xdmd2,Atild,Btild,Xp] = DMDcTest(X,gamma,B,dt);

