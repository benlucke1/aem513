function [p2p1, po2po1, M2] = normalShockRelations(M1, gamma)
%ratioPressureNormalShock - This function calculates the pressure ratios across a normal shockwave
% assumig a calorically perfect gas.

% Calculate pressure values
p2p1 = 1+ 2*(gamma)/(gamma+1).*(M1.^2 -1);
po2po1 = ((0.5*(gamma+1)*M1^2)/(1+0.5*(gamma-1)*M1^2))^(gamma/(gamma-1))*...
    (2*gamma/(1+gamma)*M1^2-(gamma-1)/(gamma+1))^(1/(1-gamma));

% Calculate mach number 
m2num = (1+((gamma-1)/2).*M1.^2);
m2den = gamma*M1.^2 - (gamma-1)/2;
M2 = sqrt(m2num/m2den);


end

