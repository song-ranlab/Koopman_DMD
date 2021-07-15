function dy = cartpend(y,qq,M,L,g,d,u)

% Sy = sin(wrapTo2Pi(y(3)));
% Cy = cos(wrapTo2Pi(y(3)));
Sy = sin(y(3));
Cy = cos(y(3));
D = qq*L*L*(M+qq*(1-Cy^2));

dy(1,1) = y(2);
dy(2,1) = (1/D)*(-qq^2*L^2*g*Cy*Sy + qq*L^2*(qq*L*y(4)^2*Sy - d*y(2))) + qq*L*L*(1/D)*u;
%  dy(3,1) = wrapTo2Pi(y(4));
dy(3,1) = y(4);
% dy(4,1) = (1/D)*((m+M)*m*g*L*Sy - m*L*Cy*(m*L*y(4)^2*Sy - d*y(2))) - m*L*Cy*(1/D)*u;
dy(4,1) = (1/D)*((qq+M)*qq*g*L*Sy - qq*L*Cy*(qq*L*y(4)^2*Sy - d*y(2))) - qq*L*Cy*(1/D)*u;