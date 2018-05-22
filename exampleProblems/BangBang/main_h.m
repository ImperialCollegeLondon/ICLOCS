% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% BangBang Control (Double Integrator Minimum Time Repositioning) Problem
%
% The problem was adapted from Example 4.11 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
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

options= settings_h(40);                  % Get options and solver settings 
[problem,guess]=BangBang;          % Fetch the problem definition

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,infoNLP1,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions
  
%%
xx=linspace(solution.T(1,1),solution.T(end,1),1000);
figure
plot(solution.T(:,1),speval(solution.Xp,1,solution.T(:,1)),'ro' )
hold on
plot(solution.T(:,1),speval(solution.Xp,2,solution.T(:,1)),'bd' )
plot(xx,speval(solution.Xp,1,xx),'r-' )
plot(xx,speval(solution.Xp,2,xx),'b-.' )
ylim([0 300])
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('States')
legend('Position [m]','Velocity [m/s]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-o' )
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u [N]')
   
    