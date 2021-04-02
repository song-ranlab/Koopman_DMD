function [A, B, nstates, augx, rXdmdc, Xdmdc] = eDMDCpoly(StateData, InputData, dt, thresh)
%Takes a full state, input vector and a time step and performs a polynomial
%regression and then DMDc on that new state. 

%% Collect and construct the snapshot 
X = StateData;
%tspan = [0.0:dt:dt*n];
[n,blah2] = size(StateData);
%%Poly Fit
usesine = 0;                                %Strictly Polynomial
polyorder = 4;                              %4th order polynomial approximation
Y = poolData(X',n,polyorder,usesine)';       %Polynomial Fitting function from Kutz et Al



%% EDMD Reconstruct

% [Phi ,omega ,lambda ,b, Xdmdc,A,B,r] = DMDcG(Y, InputData, dt, thresh);
[Phi ,omega ,lambda ,b, Xdmdc,A,B,r] = DMDcfin(Y, InputData, dt, thresh);
nstates = size(A,1);
%KO = lambda .* Phi .* inv(lambda);

%% Compute DMD mode amplitudes b
x1 = Y(:, 1);
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
Xdmdc = Phi * time_dynamics ;

rXdmdc = real(Xdmdc);
% augx = Y(1:nstates,:);
augx = Y;
end

