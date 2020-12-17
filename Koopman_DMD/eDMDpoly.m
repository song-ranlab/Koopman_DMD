function [Atilde, Xdmd,Y] = eDMDpoly(StateData, dt, thresh)
%% Collect and construct the snapshot 
X = StateData;
[m,n] = size(X);

tspan = [0.0:dt:dt*n];

%%Poly Fit
usesine = 0;                                %Strictly Polynomial
polyorder = 4;                              %4th order polynomial approximation
Y = poolData(X,m,polyorder,usesine)';       %Polynomial Fitting function from Kutz et Al
Nstates = size(Y,1);


%% EDMD Reconstruct
%% Collect and construct the snapshot 
X = Y (:,1:end -1);
Xp = Y (:,2:end);
Ups = InputData (:,1:end -1);
Omega = [Y;Ups];
%% Compute the SVD of the input space
[U,Sig ,V] = svd(Omega ,'econ');

if thresh == 0
    thresh = 1e-10;
end
    
    r = length(find(diag(Sig)>thresh));

[U, S, V] = svd(X, 'econ');
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
mm1 = size(X1, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end;
Xdmd = Phi * time_dynamics ;

