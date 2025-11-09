%%Step function
function [y] = Step(x, x0, f0, x1, f1)
if x>=x1
    y=f1;
end

if x<=x0
    y=f0;
end

if x>x0 && x<x1
    y=f0+(f1-f0)*(3*((x-x0)/(x1-x0))^2-2*((x-x0)/(x1-x0))^3);
end