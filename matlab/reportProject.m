% AEM513 - Graduate Project: 2D Diamond Airfoil Analysis
% Author: Benjamin Lucke
% Created on: 19 June 2023

%% Housekeeping
clear, clc, close all, format long
%% ################ Description & Problem Statement ##################
% This program calculates the solution to problem 5.16 in John D.
% Anderson's "Modern Compressible Flow". The problem statement is:

% Using shock-expansion theory, calculate the lift and drag (in pounds) on 
% a symmetrical diamond airfoil of semiangle ε = 15° (see Fig. 4.35) at an 
% angle of attack to the free stream of 5° when the upstream Mach number 
% and pressure are 2.0 and 2116 lb /ft2, respectively. The maximum 
% thickness of the airfoil is t = 0.5 ft. Assume a unit length of 1 ft in 
% the span direction (perpendicular to the page in Fig. 4.35).

%% ######################## Dependencies ##############################
% This program references the following User-defined functions
% betaThetaMach.m
% expansionWaves.m
% normalShockRelations.m

%% Static Pressure Calculations
%%% Constants
P_inf = 2116; % [lb/ft^2]
M_inf = 2.0; % Unitless
alfa = 5; % [deg]
epsilon = 15; % half angle of airfoil [deg]
gamma = 1.4; % Ratio of specific heats [Unitless]

%%% Leeward Side of Airfoil %%%
% We will first calculate the pressures on the compression surface and the
% expansion surface on the leeward side of the airfoil (i.e., the top
% region of the airfoil). Then, we will move to the windward side, and then
% once all pressures are calculated we may then calculate the lift and
% drag of the airfoil.

% Assuming the flow before the normal shock is isentropic,
[p1po1_top, placeholder] = isentropicRelations(M_inf);

% From the geometry, we see that the top and bottom flows are deflected at
% different angles relative to the freestream. On the top, the flow is
% deflected (15 - 5) = 10 degrees and on the bottom the flow is deflected at 
% an angle of 15+5 = 20 degrees.
turn_top = epsilon - alfa;

% The Theta-Beta-Mach relation for the weak solution gives a wave angle:
beta_top = betaThetaMach(turn_top,M_inf,1.4,1);

% Using the wave angle, we can calculate the normal component of the
% freestream velocity using Equation (4.7) and the wave angle of the 
% oblique shock:
Mn1_top = M_inf*sind(beta_top);

% Next, using the normal shock relations and the mach number of the flow
% normal to the oblique shock, we can get the properties of the
% flow across the oblique shock:
[p2p1_top,po2po1_top,Mn2_top] = normalShockRelations(Mn1_top);

% Using these properties, we can then find the static pressure on the 
% compression surface:
p_compression_top = p2p1_top*P_inf;

% We can also find the Mach number behind the shock using Equation (4.12)
% the wave angle, and the total turn angle of the flow:
M2_top = Mn2_top/sind(beta_top-turn_top);

% Next is the expansion fan. Using the Prantdl-Meyer function, Equation 
% (4.43), we can get the Prandtl-Meyer angle of the incoming flow, and
% solve for the Prandtl-Meyer angle of the downstream flow. 
v2_top = expansionWaves(M2_top);

% Since the total turn angle is 2*epsilon = 30,  
epsilon2_top = 2*epsilon;

% Equation (4.45) says that the total flow turn angle is equal to the
% Prandtl-Meyer angle of the downstream minus the PM angle of the upstream
% so theta_2 = v_3 - v_2 --> v_3 = theta_2 + v_1
v3_top = epsilon2_top + v2_top;

% Since we know the equation directly, we can iteratively solve for 
% the Mach number behind the expansion fan 
% that gives a Prandtl-Meyer angle equal to the desired value.

% Prandtl-Meyer Function (Equation (4.45))
pm_fun = @(mach) sqrt((gamma + 1)/(gamma-1)).*atan(sqrt(((gamma - 1)/...
    (gamma + 1)).*(mach.^2 - 1)))-atan(sqrt(mach.^2 - 1));

% Iteratively loop Mach numbers until error falls within .001 threshold
% from desired Prandtl-Meyer angle to arrive at downstream Mach number:
for M3_top = 1:.001:10
    v3_candidate = rad2deg(pm_fun(M3_top));
    error = abs(v3_top-v3_candidate);
    if error < .001
        break
    end
end

% The flow through an expansion fan is isentropic, so we may use the
% isentropic relations to calculate the ratio of static pressure to total
% pressure. This also means that the total pressure downstram of the
% expansion wave is identical to the total pressure upstream of the
% expansion wave.
[p3po3_top,ph] = isentropicRelations(M3_top);
p3po2_top = p3po3_top;

% Using the fact that the flow is isentropic through an expansion fan, we
% can calculate the static pressure on the expansion surface as:
p_expansion_top = p3po2_top*po2po1_top/p1po1_top*P_inf;

