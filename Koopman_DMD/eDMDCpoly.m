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

[Phi ,omega ,lambda ,b, Xdmdc,A,B] = DMDcG(Y, InputData, dt, thresh);
nstates = size(A,1);
%KO = lambda .* Phi .* inv(lambda);
rXdmdc = real(Xdmdc);
augx = Y(1:nstates,:);
end

