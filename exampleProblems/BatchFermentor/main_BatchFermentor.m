% main_h_MeshRefinement - Main script to solve the Optimal Control Problem
%
% Fed-batch fermentor optimal control problem
%
% The problem was adapted from Example 2, Section 12.4.2 of
% J.E. Cuthrell and L. T. Biegler. Simultaneous optimization and solution methods for batch reactor control profiles. Computers and Chemical Engineering, 13:49-62, 1989.
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

[problem,guess]=BatchFermentor;          % Fetch the problem definition
% options= problem.settings(200,4);
options= problem.settings(300);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figures
xx=linspace(solution.T(1,1),solution.tf,100000);
figure
hold on
plot(xx,speval(solution,'X',1,xx),'r-' )
plot(xx,speval(solution,'X',2,xx),'b-' )
plot(xx,speval(solution,'X',3,xx),'m-' )
plot(xx,speval(solution,'X',4,xx),'y-' )
plot(tv,xv(:,1),'k-.')
plot(tv,xv(:,2),'k-.')
plot(tv,xv(:,3),'k-.')
plot(tv,xv(:,4),'k-.')
xlim([0 solution.tf])
xlabel('Time [hrs]')
ylabel('States')
legend('x [g/l]','p [g/l]','s [g/l]','v [l]')
grid on

figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [hrs]')
grid on
ylabel('Control Input')
legend('u [g/hrs]')
