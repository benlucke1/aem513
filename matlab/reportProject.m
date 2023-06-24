% Title: AEM513 Graduate Project - 2D Diamond Airfoil Analysis
% Author: Benjamin Lucke
% Created on: 19 June 2023

% Revision History:
%   Repository hosted on GitHub and available on request.

%% Housekeeping
clear, clc, close all, format long

%% ~~~~~~~~~~~~~~~ Description & Problem Statement ~~~~~~~~~~~~~~~~~~~~~~
% This program calculates the solution to problem 5.16 in John D.
% Anderson's "Modern Compressible Flow". The problem statement is:

% Using shock-expansion theory, calculate the lift and drag (in pounds) on 
% a symmetrical diamond airfoil of semiangle ε = 15° (see figure in report)
% at an angle of attack to the free stream of 5° when the upstream Mach 
% number and pressure are 2.0 and 2116 lb /ft2, respectively. The maximum 
% thickness of the airfoil is t = 0.5 ft. Assume a unit length of 1 ft in 
% the span direction (perpendicular to the page in the figure).

%% ~~~~~~~~~~~~~~~~~~~~~~~~ Dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This program references the following user-defined functions:
% betaThetaMach.m
% expansionWaves.m
% normalShockRelations.m

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Constants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P_inf = 2116; % Freestream static pressure, [lb/ft^2]
M_inf = 2.0; % Freestream Mach number, Unitless
gamma = 1.4; % Ratio of specific heats [Unitless]
alfa = 5; % Airfoil angle of attack, [deg]
epsilon = 15; % Half angle of airfoil [deg]


%% Static Pressure Calculations - Top
%%% Leeward Surface of Airfoil %%%
% We will first calculate the pressures on the compression surface and the
% expansion surface on the leeward side of the airfoil (i.e., the top
% region of the airfoil). Then, we will move to the windward side, and then
% once all pressures are calculated we may then calculate the lift and
% drag of the airfoil.

% Assuming the freestream conditions are isentropic,
[p1po1_top, ~] = isentropicRelations(M_inf);

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

% Next, using the normal shock relations and the Mach number of the flow
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

% Iteratively loop Mach numbers until error falls within .01 threshold
% from desired Prandtl-Meyer angle to arrive at downstream Mach number:
for M3_top = 1:.001:100
    v3_candidate = rad2deg(pm_fun(M3_top));
    error = abs(v3_top-v3_candidate);
    if error < .01
        break
    end
end

% The flow through an expansion fan is isentropic, so we may use the
% isentropic relations to calculate the ratio of static pressure to total
% pressure after the wave. This also implies the total pressure 
% downstram of the expansion wave is identical to the total pressure 
% upstream of the expansion wave - which is the same as the total pressure
% just after the oblique shock. Hence, Po_3 = Po_2.
[p3po3_top, ~] = isentropicRelations(M3_top);
p3po2_top = p3po3_top;

% Using the fact that the flow is isentropic through an expansion fan, we
% can calculate the static pressure on the expansion surface as:
p_expansion_top = p3po2_top*po2po1_top/p1po1_top*P_inf;

%% Static Pressure Calculations - Bottom of Airfoil
%%% Windward Surface of Airfoil %%%
% We will now calculate the pressures on the compression surface and the
% expansion surface on the windward side of the airfoil (i.e., the bottom
% region of the airfoil).

% Assuming the flow before the normal shock is isentropic,
[p1po1_bottom, xxxx] = isentropicRelations(M_inf);

% From the geometry, we see that the flow is turned a total of 
% ε + alpha = 15 + 5 = 20 degrees.
turn_bottom = epsilon + alfa;

% The Theta-Beta-Mach relation for the weak solution gives a wave angle:
beta_bottom = betaThetaMach(turn_bottom,M_inf,1.4,1);

% Using the wave angle, we can calculate the normal component of the
% freestream velocity using Equation (4.7) and the wave angle of the 
% oblique shock:
Mn2_bottom = M_inf*sind(beta_bottom);

% Next, using the normal shock relations and the Mach number of the flow
% normal to the oblique shock, we can get the properties of the
% flow across the oblique shock:
[p2p1_bottom,po2po1_bottom,Mn2_bottom] = normalShockRelations(Mn2_bottom);

% Using these properties, we can then find the static pressure on the 
% compression surface:
p_compression_bottom = p2p1_bottom*P_inf;

% We can also find the Mach number behind the shock using Equation (4.12)
% the wave angle, and the total turn angle of the flow:
M2_bottom = Mn2_bottom/sind(beta_bottom-turn_bottom);

% Similar to the leeward side, use Equation (4.43) to get the Prandtl-Meyer
% angle for the Mach number upstream of the fan. 
v2_bottom = expansionWaves(M2_bottom);

% Again, the total turn angle is 2*epsilon = 30,  
epsilon2_bottom = 2*epsilon;

% Equation (4.45) says that the total flow turn angle is equal to the
% Prandtl-Meyer angle of the downstream minus the PM angle of the upstream
% so theta_2 = v_3 - v_2 --> v_3 = theta_2 + v_1
v3_bottom = epsilon2_bottom + v2_bottom;

% Identical to above, we can iteratively solve for the Mach number that 
% yields the correct Prandtl-Meyer angle for the bottom expansion fan:
for M3_bottom = 1:.001:100
    v3_candidate = rad2deg(pm_fun(M3_bottom));
    error = abs(v3_bottom-v3_candidate);
    if error < .01
        break
    end
end

% The flow through an expansion fan is isentropic, so we may use the
% isentropic relations to calculate the ratio of static pressure to total
% pressure after the wave. This also implies the total pressure 
% downstram of the expansion wave is identical to the total pressure 
% upstream of the expansion wave - which is the same as the total pressure
% just after the oblique shock. Hence, Po_3 = Po_2.
[p3po3_bottom, ~] = isentropicRelations(M3_bottom);
p3po2_bottom = p3po3_bottom;

% Using the fact that the flow is isentropic through an expansion fan, we
% can calculate the static pressure on the expansion surface as:
p_expansion_bottom = p3po2_bottom*po2po1_bottom/p1po1_bottom*P_inf;

%% Calculation of Lift and Drag

% The pressure forces act as a distributed load on the surfaces of the
% airfoil in 2 dimensions. To calculate the magnitude of the forces (i.e.
% lift and drag), we must multiply the area of the airfoil's surface by 
% the the Pressure forces on each surface, and then resolve the components
% in a coordinate system.

% The maximum thickness of the airfoil was given as 0.5 feet, and the
% spanwise length of the airfoil was given as 1 foot.
t = 0.5; % Airfoil thickness, [ft]
w = 1; % Airfoil span, [ft]

% From the geometry, we see that the half chord length of the airfoil is
% give by the quantity c/2 and the hypotenuse, l, is not given. Using basic
% trigonometry, the hypotenuse length of the triangle is calculated as half
% of the thickness of the airfoil divided by the sin of the semispan angle,
% epsilon.
l = (t/2)/sind(epsilon);

% From the figure, we see the area of the surfaces being acted on by the
% pressure distributions is a square, and thus the area of each compression and
% expansion surface is given as l*w
area_surfs = l*w;

% Thus, the magnitude of the normal force vectors are:
F_compression_top_mag = p_compression_top*area_surfs;
F_expansion_top_mag = p_expansion_top*area_surfs;
F_compression_bottom_mag = p_compression_bottom*area_surfs;
F_expansion_bottom_mag = p_expansion_bottom*area_surfs;

% The coordinate system of choice in this analysis is the following:
%                  
%                       ^ L+
%                       |
%                       |
%                       |
%                       ----------> D+
%

% By inspection of the geometry, the components of the normal force vectors
% for drag and lift are found respectively as D = F*sin() and L = F*cos():
F_compression_top = [F_compression_top_mag*sind(epsilon-alfa); 
                    -F_compression_top_mag*cosd(epsilon-alfa)];

F_expansion_top = [-F_expansion_top_mag*sind(epsilon+alfa); 
                    -F_expansion_top_mag*cosd(epsilon+alfa)];

F_compression_bottom = [F_compression_bottom_mag*sind(epsilon+alfa); 
                    F_compression_bottom_mag*cosd(epsilon+alfa)];

F_expansion_bottom = [-F_expansion_bottom_mag*sind(alfa); 
                    F_expansion_bottom_mag*cosd(alfa)];

L = F_compression_top(2) + F_expansion_top(2) + F_compression_bottom(2) + F_expansion_bottom(2);
D = F_compression_top(1) + F_expansion_top(1) + F_compression_bottom(1) + F_expansion_bottom(1);


%% GENERATING TABLES

T = table([F_compression_top; norm(F_compression_top)],...
    [F_expansion_top; norm(F_expansion_top)] ...
    ,[F_compression_bottom; norm(F_compression_bottom)],...
    [F_expansion_bottom; norm(F_expansion_bottom)],...
    [D; L; norm(F_compression_top)+norm(F_expansion_top)+norm(F_expansion_bottom)+norm(F_compression_bottom)],...
    'RowNames',{'Drag (lb)','Lift (lb)','Magnitude (lb)'},...
    'VariableNames',{'Top Compression Surface','Top Expansion Surface',...
    'Bottom Compression Surface','Bottom Expansion Surface','Airfoil'})
