% Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Time optimal control of a two-link robot arm
%
% The problem was adapted from Example 2, Section 12.4.2 of
% Rein Luus. Iterative Dynamic Programming. Chapman & Hall/CRC Monographs and Surveys in Pure and Applied Mathematics, 2000.
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

[problem,guess]=TwoLinkRobotArm;          % Fetch the problem definition
options= problem.settings(20);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figure

xx=linspace(solution.T(1,1),solution.T(end,1),1000);
figure
plot(xx,speval(solution,'X',1,xx),'r-' )
hold on
plot(xx,speval(solution,'X',2,xx),'b-' )
plot(xx,speval(solution,'X',3,xx),'m-' )
plot(xx,speval(solution,'X',4,xx),'g-' )
plot(tv,xv(:,1),'k-.')
plot(tv,xv(:,2),'k-.')
plot(tv,xv(:,3),'k-.')
plot(tv,xv(:,4),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('States')
legend('\omega_{\phi} [rad/s]','\omega_{\psi} [rad/s]','\chi [rad]','\phi [rad]')
grid on

figure
plot(xx,speval(solution,'U',1,xx),'b-' )
hold on
plot(xx,speval(solution,'U',2,xx),'r-' )
plot(tv,uv(:,1),'k-.')
plot(tv,uv(:,2),'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u1','u2')
