function [A, B, nstates, Y, rXdmdc, Xdmdc] = eDMDCFFT(StateData, InputData, dt, order, thresh)
%Takes a full state, input vector and a time step and performs a polynomial
%regression and then DMDc on that new state. 

% A: Learned A matrix
% B: Learned B matrix
% nstates: Number of states in learned matrix
% augx: lifted states truncated to DMD truncation 
% rXdmdc: Just real values of Xdmdc reconstruction - Pretty useless
% Xdmdc: Complex of Xdmdc reconstruction - not as useless

%% Collect and construct the snapshot 
X = StateData;
%tspan = [0.0:dt:dt*n];
[n,blah2] = size(StateData);
%%Poly Fit
Y = fft(StateData,order,1);



%% EDMD Reconstruct

[Phi ,omega ,lambda ,b, Xdmdc,A,B] = DMDcG(Y, InputData, dt, thresh);
nstates = size(A,1);
%KO = lambda .* Phi .* inv(lambda);
rXdmdc = real(Xdmdc);
%augx = Y(1:nstates,:);
end

