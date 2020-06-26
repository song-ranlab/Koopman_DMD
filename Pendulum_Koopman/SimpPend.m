function [dy] = SimpPend(y,u)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dy =([y(2); -sin(y(1))+u]);
end

