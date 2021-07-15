function [Phi ,omega ,lambda ,b,Xdmd,approxA,approxB,r] = DMDcfin(LiftData, InputData, dt, thresh)
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

%% COmpute the SVD of the input space
[U,Sig ,V] = svd(Omega ,'econ');
%p = length(find(diag(Sig)>thresh));        %Input space truncation 
if thresh == 0
    p = size(Omega,1);
else
    if thresh < 0
        p = length(find(diag(Sig)>-thresh));
   
    elseif thresh >= 0 && thresh <= size(Omega,1)
        p = floor(thresh);
        if p == 0
            p = size(Omega,1);
        end
    else
        
        p = size(Omega,1)
        %p = length(find(diag(Sig)>thresh));
    end
   
end
Utilda = U(:,1:p);
Sigtilda = Sig(1:p ,1:p);
Vtilda = V(:,1:p);

%% Compute the SVD of the output space X'
[U,Sig ,V] = svd(Xp,'econ');

if thresh == 0
    r = size(Xp,1);
else
    if thresh <0
        r = length(find(diag(Sig)>-thresh));
    elseif thresh >=0 && thresh <= size(Omega,1)
        r = floor(thresh);
        if r == 0
            r = size(Omega,1);
        end
    else
        r = size(Omega,1)
        %p = length(find(diag(Sig)>thresh));
    end
   
end
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
% Atilde = U_r' * X2 * V_r / S_r; % low-rank dynamics

% checktilda = inv(Sigtilda);

% approxA = (Xp - [[0;.2;0;.1];zeros((size(Xp,1)-4),1)]*Ups)*Vbar*inv(Sighat)*Uhat'
% approxA = (Xp)*Vtilda*inv(Sigtilda)*U1';
approxB = Uhat'*(Xp)*Vtilda*inv(Sigtilda)*U2';
% AB = Xp*Vtilda*Sigtilda^(-1)*Utilda';
% approxA = AB(1:Nstates,1:Nstates);
% approxB = AB(1:Nstates,end); 



% approxA = (Xp - approxB*Ups)*V*inv(Sig)*U';
%% Preform the eigenval decomp of Atilda
[W,D] = eig(approxA);

%Compte the dynamic modes of the operator A
Phi = Xp * Vtilda * inv(Sigtilda) * U1'*Uhat * W;
% Phi = Xp * Vtil * inv(Sigtil) * U_2'*Uhat * W;
% Finish to output
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
% x1 = Omega(:, 1);
x1 = X(:,1);
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
% AB = Xp*Vtilda*Sigtilda^(-1)*Utilda';
% approxA = AB(1:Nstates,1:Nstates);
% approxB = AB(1:Nstates,end); 
end

