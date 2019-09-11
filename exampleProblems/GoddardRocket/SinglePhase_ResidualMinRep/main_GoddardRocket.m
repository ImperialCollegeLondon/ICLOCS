% main_GoddardRocket - Main script to solve the Optimal Control Problem
%
% Goddard Rocket Problem (Single-phase Formulation)
%
% The problem was adapted from Example 4.9 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


%--------------------------------------------------------




%% Solve the problem
clear all;close all;format compact;
options= settings_GoddardRocket(100);                  % Get options and solver settings 
[problem,guess]=GoddardRocket;          % Fetch the problem definition
[solution,MRHistory]=solveMyProblem( problem,guess,options);
genSolutionPlots(options, solution);

%% Generate Figures
xx=linspace(solution.T(1,1),solution.tf,100000);

figure(100)
hold on
plot(xx,speval(solution,'X',1,xx))
xlabel('Time [s]')
ylabel('Altitude [km]')
grid on

figure(101)
hold on
plot(xx,speval(solution,'X',2,xx))
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on

figure(102)
hold on
plot(xx,speval(solution,'X',3,xx))
xlabel('Time [s]')
ylabel('Mass [kg]')
grid on

figure(103)
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
plot(xx,speval(solution,'U',1,xx))    
ylim([problem.inputs.ul problem.inputs.uu])
xlabel('Time [s]')
grid on
ylabel('Control Input')