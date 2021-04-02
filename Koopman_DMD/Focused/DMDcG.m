function [Phi ,omega ,lambda ,b,Xdmd,approxA,approxB,r] = DMDcG(LiftData, InputData, dt, thresh)
%DMDC
% function [Phi,omega ,lambda ,b,Xdmd ] = DMD(X1,X2,r,dt)
% Computes the Dynamic Mode Decomposition of State Data with Control input
%
% INPUTS:
% State Input - All measurements from k = 0 to m timesteps
% Input Data - All control inputs from k = 0 to m-1 (last input and first
% input should be zero
% dt = time step advancing X to X’
% thresh - Threshold for DMD ranks, currently set to same input value, if no value is given set to 1e-10 
%
% OUTPUTS:
% Phi - the DMD modes
% omega - the continuous-time DMD eigenvalues
% lambda - the discrete-time DMD eigenvalues
% b - a vector of magnitudes of modes Phi
% approxA - 
% Xdmd - the data matrix reconstructed by Phi, omega , b
%% Collect and construct the snapshot 
X = LiftData (:,1:end -1);
Xp = LiftData (:,2:end);
Ups = InputData (:,1:end -1);
Omega = [X;Ups];

 if thresh == 0
     thresh = 1e-10;
 end

%% COmpute the SVD of the input space
[U,Sig ,V] = svd(Omega ,'econ');
p = length(find(diag(Sig)>thresh));        %Input space truncation 
% p = size(Omega,1);
Utilda = U(:,1:p);
Sigtilda = Sig(1:p ,1:p);
Vtilda = V(:,1:p);

%% Compute the SVD of the output space X'
[U,Sig ,V] = svd(Xp,'econ');
%if thresh == 0
    r = length(find (diag(Sig)>thresh));        % Output space truncation
%else
% r=p-1;
    
Uhat = U(:,1:r);
Sighat = Sig(1:r,1:r);
Vbar = V(:,1:r);

%Compute A & B
n = size(X,1);
q = size(Ups ,1);
U1 = Utilda (1:n,:);
U2 = Utilda (n+q:n+q,:);
Nstates = size(LiftData,1);
approxA = Uhat'*(Xp)*Vtilda*inv(Sigtilda)*U1'*Uhat;
approxB = Uhat'*(Xp)*Vtilda*inv(Sigtilda)*U2';
% AB = Xp*Vtilda*Sigtilda^(-1)*Utilda';
% approxA = AB(1:Nstates,1:Nstates);
% approxB = AB(1:Nstates,end);



% approxA = (Xp - approxB*Ups)*V*inv(Sig)*U';
%% Preform the eigenval decomp of Atilda
[W,D] = eig(approxA);

%Compte the dynamic modes of the operator A
Phi = Xp * Vtilda * inv(Sigtilda) * U1'*Uhat * W;

%Finish to output
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
x1 = X(:, 1);
%x1 = X(:,:);
b = Phi\x1;%double check ifft here
%put in euler prop
%
mm1 = size(X, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1,
    time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
Xdmd = Phi * time_dynamics ;
      
end

