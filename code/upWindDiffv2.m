function [difference] = upWindDiffv2(f,u,dx)
    difff = (f(2:end) - f(1:end-1))/dx; 
    backdiff = [(f(1) - f(end))/dx, difff];    % backwards difference
    fordiff = [difff,(f(1) - f(end))/dx];      % forwards difference
    posVel = u > 0;                            % locations of positive velocity
    difference = fordiff.*~posVel + backdiff.*posVel;   
end

