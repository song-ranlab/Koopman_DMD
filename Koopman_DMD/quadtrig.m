function [dy] = quadtrig(inputstate, controlinput, characteristics)
x = inputstate;
u = controlinput;
g = -9.81;
I_x = characteristics(1);
I_y = characteristics(2);
I_z = characteristics(3);
m = characteristics(4);

g1x = (-1/m * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))));
g2x = (-1/m * (sin(x(7))*cos(x(9)) - cos(x(7))*sin(x(9))*sin(x(8))));
g3x = (-1/m * (cos(x(7))*cos(x(8))));

dy(1,1)= x(6) * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) + x(5) * (sin(x(7))*cos(x(9))*sin(x(8)) - sin(x(9))*cos(x(7))) + x(4) * (cos(x(8))*cos(x(9)));
dy(2,1)= x(5) * (cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(9))*sin(x(8))) + x(6) * (sin(x(7))*cos(x(9))*sin(x(8)) - cos(x(9))*sin(x(7))) + x(4) * (cos(x(8))*sin(x(9)));
dy(3,1)= x(6) * (cos(x(7))*cos(x(8))) - x(4) * sin(x(8)) + x(5) * (cos(x(8))*sin(x(7))); 
dy(4,1)= -u(1)/m * (sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) + g1x*u(1);
dy(5,1)= -u(1)/m * (cos(x(7))*sin(x(9))*sin(x(8)) - cos(x(9))*sin(x(7))) + g2x*u(1);
dy(6,1)= g - u(1)/m * (cos(x(8))*cos(x(9))) + g3x*u(1);
dy(7,1)= x(10) + x(12) * (cos(x(7))*tan(x(8))) + x(11) * (sin(x(7))*tan(x(8)));
dy(8,1)= x(11)*cos(x(7)) - x(12)*sin(x(7));
dy(9,1)= x(12)*(cos(x(7))/cos(x(8))) + x(11)*(sin(x(7))/cos(x(8)));
dy(10,1)= ((I_y-I_z)/I_x) * x(11)*x(12) + u(2)/I_x;
dy(11,1)= ((I_z-I_x)/I_y) * x(10)*x(12) + u(3)/I_y;
dy(12,1)= ((I_x-I_y)/I_z) * x(10)*x(11) + u(4)/I_z;
                 
end

