function [xhist,uhist,thist,rhist] = HBMPC(syshandle, Duration,Ts,N,Nu,x0, Q, R, Ru, xref, pest,uval)
global order
%% Definitions
%Syshandle: Target System
%Duration: Duration of maneuver
%Ts: Time Step
%N: Length of prediction Horizon
%Nu: Control Horizon
%ObjectiveFCN: Function call for system we intend to optimize
%ConstraintFCn: Constraints we wish to impose on the system

global Nvar order
% Prepare variables
Nvar = length(xref);
Nt = (Duration/Ts)+1;
%uopt0    = firstu;
uopt0    = 0;
xhat     = x0;
uopt     = uopt0.*ones(Nu,1);
xHistory = zeros(Nvar,Nt); 
xHistory(:,1) = xhat;
uHistory = zeros(1,Nt);    
uHistory(1)   = uopt(1);
tHistory = zeros(1,Nt);    
tHistory(1)   = 0;
rHistory = zeros(Nvar,Nt);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
% Start simulation
fprintf('Chinchilla done snacking. Run wheel engaged...\n')
tic
for ct = 1:Nt-1
    % NMPC with full-state feedback
    COSTFUN = @(u) ObjectiveFCN(u,xhat,N,Nu,xref,uHistory(:,ct),pest,Q,R,Ru);
    uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],[],options);
    
    % Integrate system
    xhat = myrk4u(syshandle,xhat,uopt(1),Ts/10,10,[],[]);  
%     [~,xhat] = ode45(syshandle,[0:Ts:Duration],uopt(1))
    xhist(:,ct+1) = xhat;
    uhist(:,ct+1) = uopt(1);
    thist(:,ct+1) = ct*Ts;
    rhist(:,ct+1) = xref;
    
    if mod(ct,10) == 0
        disp(['Chinchilla Endurance: ',num2str(100*ct/(Duration/Ts)),'% tired'])
    end
end



%% Objective Function Definition

function J = ObjectiveFCN(uopt,x,N,Nu,xref,u0,p,Q,R,Ru)

global Nvar order
u = uopt;
%checkers=poolData(x',Nvar,p.polyorder,p.usesine)'
%% Integrate system
% [xk,~] = lsim(p.sys,[u' 0],[0:N].*p.sys.Ts,poolData(x',Nvar,p.polyorder,p.usesine)');
[xk,~] = lsim(p.sys,[u' 0],[0:N].*p.sys.Ts,fft(x,order,1));
 xk = xk';
% xk = xk(1:p.nvar,:);

%% Cost Calculation
% Set initial plant states, controller output and cost.
uk = u(1);
J = 0;
% Loop through each prediction step.
for ct=1:N
    % Obtain plant state at next prediction step.
    xk1 = xk(:,ct);

    % accumulate state tracking cost from x(k+1) to x(k+N).
%     J = J + (p.Hfun(xk1)-p.Hfun(xref))'*Q*(p.Hfun(xk1)-p.Hfun(xref));
%     J = J + (xk1-poolData(xref',Nvar,p.polyorder,p.usesine)')'*Q*(xk1-poolData(xref',Nvar,p.polyorder,p.usesine)');
        J = J + (xk1-fft(xref,order,1))'*Q*(xk1-fft(xref,order,1));
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J = J + ((uk-u0)'*R*(uk-u0)) + uk'*Ru*uk;
    else
        J = J + ((uk-u(ct-1))'*R*(uk-u(ct-1))) + uk'*Ru*uk;
    end
    % Update uk for the next prediction step.
    if ct<N
        uk = u(ct+1);
    end
end


%% 4th order Rune-Kutta

function x = rk4u(f,x,u,h,n,t,p)

% RK4U   Runge-Kutta scheme of order 4 for control system
%   rk4u(v,X,U,h,n) 
%n = number of steps taken
%h = stepsize
%x = initial condition
%t & p are parameters for function calls of f(t,x,u,p) standard
%   f(X,U) dynamical system function handle

for i = 1:n
    k1 = f(t,x,u,p); 
    k2 = f(t,x + h/2*k1,u,p);
    k3 = f(t,x + h/2*k2,u,p);
    k4 = f(t,x + h*k3,u,p);
    x = x + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end

