% Main script to solve the Optimal Control Problem 
%
% Aly-Chan problem
%
% The problem was originally presented by: 
% G.M. Aly and W.C. Chan. 1973.  Application of a modified quasilinearization technique to totally singular optimalcontrol problems.Internat. J. Control17, 4 (1973), 809?815
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

[problem,guess]=AlyChan;          % Fetch the problem definition
options= problem.settings(100);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figure

xx=linspace(solution.T(1,1),solution.T(end,1),5000);
figure
hold on
plot(xx,speval(solution,'X',1,xx),'b-' )
plot(xx,speval(solution,'X',2,xx),'r-' )
plot(xx,speval(solution,'X',3,xx),'g-' )
plot(tv,xv(:,1),'k-.' )
plot(tv,xv(:,2),'k-.' )
plot(tv,xv(:,3),'k-.' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('States')
legend('x1 [-]','x2 [-]','x3 [-]')
grid on

figure
plot(xx,speval(solution,'U',1,xx),'b-' )
hold on
plot(tv,uv(:,1),'k-.' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlim([0 solution.tf])
ylim([-1 1])
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u [N]')
   
