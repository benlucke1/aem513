function [p2p1, po2po1, M2] = normalShockRelations(M1, gamma)
%NORMALSHOCKRELATIONS This function calculates the pressure ratios across a normal shockwave
% assuming a calorically perfect gas.
%
% Inputs
%  M1 - Upstream Mach Number [Real]
%  gamma - Ratio of specific heats [Real]
%
% Outputs
%  p2p1 - Static pressure ratio across normal shock
%  po2po1 - Stagnation (total) pressure ratio across normal shock
%  M2 - Downstream Mach Number across normal shock

if nargin < 2
    gamma = 1.4;
end
% Calculate pressure values
p2p1 = 1+ 2*(gamma)/(gamma+1).*(M1.^2 -1);
po2po1 = ((0.5*(gamma+1)*M1^2)/(1+0.5*(gamma-1)*M1^2))^(gamma/(gamma-1))*...
    (2*gamma/(1+gamma)*M1^2-(gamma-1)/(gamma+1))^(1/(1-gamma));

% Calculate mach number 
m2num = (1+((gamma-1)/2).*M1.^2);
m2den = gamma*M1.^2 - (gamma-1)/2;
M2 = sqrt(m2num/m2den);


end

