% MAIN - Main script to solve the Optimal Control Problem
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all;format compact;


[problem,guess]=myProblem;  		% Fetch the problem definition

% options= settings_h(40);                  % Get options and solver settings 
% options= settings_hp(5,4);                  % Get options and solver settings 
options= settings_hp([4 5 3],[-1 0.3 0.4 1]);                  % Get options and solver settings 

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions

xx=linspace(solution.T(1,1),solution.tf,1000);
% x1=speval(solution.Xp,1,xx);
x1=speval(solution.Xp,1,solution.TSeg_Bar,xx);

