% MAIN - Main script to solve the Optimal Control Problem
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
clear all;close all;format compact;

[problem,guess]=myProblem;          % Fetch the problem definition
options= problem.settings(Nps);                  % Get options and solver settings (h method)
% options= problem.settings(Nps,Npd);            % Get options and solver settings (hp method)
% options= problem.settings([-1,...,1],[Npd_1,...Npd_N]);            % Get options and solver settings (hp method)
[solution,MRHistory]=solveMyProblem( problem,guess,options); % Format for NLP solver, solve and post-process the NLP


%% Optional functions

% Built-in figure generation (check settings for the type of plots to be generated)
genSolutionPlots(options, solution); 

% Simulate the obtained solution (open-loop) with ODE integration 
% [ tv, xv, uv ] = simulateSolution( problem, solution);
% [ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113');
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 'ode113', 0.01 );

% Obtainting the interpolated solution from polynomials
tt=linspace(solution.t0,solution.tf,1000);
x1=speval(solution,'X',1,tt);
x2=speval(solution,'X',2,tt);
...
u1=speval(solution,'U',1,tt);
u2=speval(solution,'U',2,tt);
...