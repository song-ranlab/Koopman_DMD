function [y] = pendulumeom(x,u)
y(1,1) = x(2); 
y(2,1) = -sin(x(1))+u;
end

