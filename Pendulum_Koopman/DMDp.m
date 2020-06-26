function [Phi ,omega ,lambda ,b,Xdmd] = DMDp(X1,X2,r,dt,lspan,tend)
% function [Phi,omega ,lambda ,b,Xdmd ] = DMD(X1,X2,r,dt)
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS:
% X1 = X, data matrix
% X2 = X’, shifted data matrix
% Columns of X1 and X2 are state snapshots
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X’)
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous -time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstructed by Phi, omega , b
%% DMD
[U, S, V] = svd(X1, 'econ');
r = min(r, size(U,2));
U_r = U(:, 1:r); % truncate to rank -r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * X2 * V_r / S_r; % low-rank dynamics
[W_r , D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
x1 = X1(:, 1);
b = Phi\x1;
%% DMD reconstruction
%mm1 = size(X1, 2); % mm1 = m - 1
mm1 = lspan;
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
% Xdmd = Phi * time_dynamics ;
%% DMD Prediction
mm2 = tend;
time_dynamics2 = zeros(r, (mm2-lspan));
t = (mm1:mm2 -1)*dt; % time vector
for iter = mm1+1:mm2 ,
time_dynamics2 (:,(iter-mm1)) = (b.*exp(omega*t(iter)));
end;
% Xdmd = Phi * time_dynamics2 ;
time_dynamics3=zeros(r,mm2);
for i = 1:mm1
    time_dynamics3(:,i)=time_dynamics(:,i);
end
for i = 1:(mm2-mm1)
    time_dynamics3(:,i+mm1)=time_dynamics2(:,i);
end
Xdmd = Phi * time_dynamics3 ;

