function [p_p0, t_t0] = isentropicRelations(M, gamma)
%ISENTROPICRELATIONS Calculates the ratio of static to stagnation pressure
%and temperatures of a flow given a Mach number assuming isentropic flow
%and calorically perfect gas.
%
% Inputs
%  M1 - Mach Number [Real]
%  gamma - Ratio of specific heats [Real] (default = 1.4)
%
% Outputs
%  p_p0 - Static to stagnation (total) pressure ratio (P/P_0)
%  t_p0 - Static to stagnation (total) temperature ratio (T/T_0)

if nargin < 2
    gamma = 1.4;
end
g = gamma;

p_p0 = (1+((g-1)/2).*M.^2).^((-1*g)./(g-1));
t_t0 = (1+((g-1)/2).*M.^2).^(-1);
end


%% SOURCE OF EQUATIONS:
% https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
