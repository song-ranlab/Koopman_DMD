function [y] = GenLinData(A,B,x0,u,duration )
% Solve a linear ss using an initial condition and control inputs
m = size(duration);
[n,l] = size(A);
% y = nan(n,m(2));
y(:,1) = x0;
xkp1 = x0;

for i = 1:m(2)-1
   y(:,i+1) =  A*xkp1 + B*u(i);
   xkp1 = y(:,i+1);
end

end

