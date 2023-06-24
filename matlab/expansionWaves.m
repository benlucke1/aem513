function v = expansionWaves(M1, gamma)

%NORMALSHOCKRELATIONS This function calculates the Prandtl-Meyer angle
%across an oblique shock assuming a calorically perfect gas.
%
% Inputs
%----------
%  M - Local Mach Number [Real]
%  gamma - Ratio of specific heats [Real] (default = 1.4)
%
% Outputs
%----------
%  v - Prandtl-Meyer angle

%% Error Handling 
if nargin < 2
    gamma = 1.4;
end
if nargin < 1
    error('Error. Please input the Mach number.')
end

%% MAIN
g = gamma;
m = M1;

t1 = sqrt((gamma + 1)/(gamma-1));
t2 = atan(sqrt(((gamma - 1)/(gamma + 1)).*(m^2 - 1)));
t3 = atan(sqrt(m^2 - 1));

v = rad2deg(t1*t2-t3);
