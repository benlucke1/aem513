function [beta, theta, mach] = betaThetaMach(theta, mach, gamma, delta)
%BETATHETAMACH This function calculates the wave angle, beta, of an oblique
% shock given the deflection angle, theta (in degrees) and the upstream
% Mach number
%   Inputs:
%       theta - deflection angle in degrees
%       mach - upstream mach number of flow
%       gamma - specific heat ratio, typically 1.4
%       delta - [0,1] - 0 is strong shock solution, 1 is weak shock
%       solution
%
%   Outputs:
%       beta - wave angle of oblique shock in degrees
%       theta - input deflection angle in degrees
%       mach - input mach number
t = theta;
m = mach;
del = delta;
g = gamma;


lam = sqrt(((m.^2 - 1).^2)-3*(1+(g-1)/2.*m.^2).*(1+(g+1)/2.*m.^2).*(tand(theta).^2));
chi_n = ((m.^2 - 1).^3) - 9*(1+(g-1)/2.*m.^2)*(1+((g-1)/2).*m.^2+((g+1)/4).*m.^4)*(tand(theta).^2);
chi_d = lam.^3;
chi = chi_n./chi_d;
tanb_num = (m.^2 - 1 + 2*lam*cos((acos(chi)+ 4*pi*delta)/3));
tanb_den = 3*(1+((g-1)/2).*m.^2)*tand(theta);
tanb = tanb_num./tanb_den;
b = rad2deg(atan(tanb));


beta = b;
theta = theta;
mach = mach;
end

