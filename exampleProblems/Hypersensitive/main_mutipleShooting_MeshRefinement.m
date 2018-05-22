% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Hypersensitive problem
%
% The problem was originally presented by: 
% A.V. Rao, and K.D. Mease, "Eigenvector approximate dichotomic basis method for solving hyper?sensitive optimal control problems", Optimal Control Applications and Methods, 21(1), pp.1-19, 2000
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

clear all;close all;format compact;
global sol;  
sol=[];                             % Initialize solution structure

options= settings_multipleShooting(40);                  % Get options and solver settings 
[problem,guess]=Hypersensitive;          % Fetch the problem definition
[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,infoNLP1,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);         % Output solutions
    
%%
xx=linspace(solution.T(1,1),solution.tf,1000*max(1000,length(solution.X)));

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'ro' )
hold on
plot(xx,speval(solution.Xp,1,xx),'r-' )
xlabel('Time [s]')
ylabel('x(t)')
grid on


figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-o' )
xlabel('Time [s]')
grid on
ylabel('u(t)')

