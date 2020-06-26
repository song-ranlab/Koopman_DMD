function [Phi ,omega ,lambda ,b,Xdmd] = DMDc(X1,X2,U1,U2,OMEGA,r,p,B,dt)
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
% [Uu, Su, Vu] = svd(U1, 'econ');
[U3, S2, V2] = svd(X2, 'econ');
[Uh, Sh, Vh] = svd(OMEGA, 'econ');
r = min(r, size(U,2));
U_r = U(:, 1:r); % truncate to rank -r
U_r2 = U3(:, 1:r); % truncate to rank -r
U_h = Uh(1:2, 1:r); % truncate to rank -r
S_r = S(1:r, 1:r);
S_t = Sh(1:p, 1:p);
V_r = V(:, 1:r);
V_t = Vh(:, 1:p);
Vc = V2(:,1:r);
Sig  = S2(1:r,1:r);
Ue = U_h' * U_r2;
Pb = V_t / S_t;
Atilde =U_r2' * X1 * Pb * Ue ; % low-rank dynamics
% Atilde2 = (X2-B*U2)*V*inv(Sig)*U_r2')
% Atilde =U_r2' * X2 * V_t / S_t * U_h' * U_r2; % low-rank dynamics
[W_r , D] = eig(Atilde);
Phi = X1 * V_t / S_r * U_h' * U_r2 * W_r; % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
x1 = X1(:, 1);
b = Phi\x1;
%% DMD reconstruction
mm1 = size(X1, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
Xdmd = Phi * time_dynamics ;

%% Other Algorithm from demo
% Vc = V2(:,1:r);
% Sig  = S2(1:r,1:r);
A = ((X2-B*U2)*V_r*inv(Sig)*U_r2');
[W_r , D] = eig(A);
Phi = X1 * V_t / S_r * U_r2 * W_r; % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
x1 = X1(:, 1);
b = Phi\x1;
%% DMD reconstruction
mm1 = size(X1, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
Xdmd = Phi * time_dynamics ; 
 
 