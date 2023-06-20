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


%% Solution
% We will first calculate the pressures on the compression surface and the
% expansion surface on the leeward side of the airfoil (i.e., the top
% region of the airfoil)
alfa = 5;
Minf = 2;
epsilon = 15;

[beta_1_top, theta_1_top, m1] = betaThetaMach(epsilon-alfa,Minf,1.4,1);
mn1 = m1*sind(beta_1_top);
[p2p1, tp2tp1, M2] = ratioPressureNormalShock(mn1,1.4);
