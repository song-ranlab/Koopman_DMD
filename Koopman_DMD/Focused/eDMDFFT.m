function [Atilde, Xdmd,Y] = eDMDFFT(StateData, dt, order, thresh)
%% Collect and construct the snapshot 
%X = StateData;
[m,n] = size(StateData);

tspan = [0.0:dt:dt*n];

% Y = fft2(StateData);


% Y = fft(StateData,order,1); 
Y = fft(StateData',n); 
Y=Y';

%Y = Y';
% PSD = Y.*conj(Y)/length(tspan);
% fz = 1/(dt*n)*(0:n);
% L = 1:floor(n/2);

%figure,plot(fz(L),PSD(L)),hold on
% % 
% indices = PSD>10;
% PSDclean = PSD.*indices;
% Y = indices.*Y;
% figure,plot(fz(L),PSDclean(L)),hold on
% Y2 = ifft(Y);
% figure,plot(Y2(1,1:n)),hold on, plot(Y2(2,1:n))
% Fs = 1/dt;
% T = dt;
% L = length(tspan);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot
% Y = P2';
% Y = Ynew(:,:,1);
% Y = fft(StateData',order,2); %THIS IS WRONG DONT DO IT! JIM, DONT DO ITT!
%% Iterative


% for j=1:size(StateData,1);  % Compute row-wise FFT
%     Cshift(j,:) = fftshift(fft(StateData(j,:)));
%     C(j,:) = (fft(StateData(j,:)));
% end
% 
% for j=1:size(C,2);  % Compute column-wise FFT
%     Y(:,j) = fft(C(:zo,j));
% end
% % Fourier Fit
% 
% set x and y
% xo1 = StateData(1,1:end-1);
% yo1 = StateData(1,2:end);
% [Ao, Bo,Y1] = Fseries(xo1,yo1,4,false);
% 
% xo2 = StateData(2,1:end-1);
% yo2 = StateData(2,2:end);
% [Ao2, Bo2,Y2] = Fseries(xo2,yo2,4,false);
% 
% xo = StateData(:,1:end-1);
% yo = StateData(:,2:end);
% [Ao, Bo,Y] = Fseries(xo',yo',4,false);
% 
% Y1 = Fseriesval(Ao,Bo,tspan);
% Y2 = Fseriesval(Ao2,Bo2,tspan);

%% Collect and construct the snapshot 
X = Y(:,1:end -1);
Xp = Y (:,2:end);
%Ups = InputData (:,1:end -1);
%Omega = [Y;Ups];
%% Compute the SVD of the input space
[U,Sig ,V] = svd(X ,'econ');

if thresh == 0
   r = length(find(diag(Sig)>thresh));
else
    r = thresh;
end
    

[U, S, V] = svd(X, 'econ');
%r = min(r, size(U,2));
U_r = U(:, 1:r); % truncate to rank -r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * Xp * V_r / S_r; % low-rank dynamics
[W_r , D] = eig(Atilde);
Phi = Xp * V_r / S_r * W_r; % DMD modes
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
%x1 = poolData(X(1, :),m,polyorder,usesine)';
x1 =  X(:, 1);
% b = Phi\ifft(x1,16,1)
b = Phi\x1;
%% DMD reconstruction
mm1 = size(X, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:mm1 ,
time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd2 = Phi * time_dynamics ;
Xdmd = Xdmd2';
% fhat = Phi * time_dynamics;
% psd = fhat.*conj(fhat)/n;
% freq = 1/(dt*n)*(0:n);
% L=1:floor(n/2);
% 
% indicies = psd>100;
% Xdmd = indicies.*fhat;
% pause
