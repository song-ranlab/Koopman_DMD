function dy = mycartpend(y,m,M,L,g,d,u)

Sy = sin(y(3));
Cy = cos(y(3));
D = m*L*L*(M+m*(1-Cy^2));

dy(1,1) = y(2);
dy(2,1) = (u + m*Sy*(L*y(4)^2-g*Cy))/(M+m*Cy^2);
dy(3,1) = y(4);
dy(4,1) = (-u*Cy - m*L*y(4)^2*Cy*Sy+(M+m)*g*Sy)/(L*(M+m*Sy^2));