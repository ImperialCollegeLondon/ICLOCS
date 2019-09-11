% Main script to solve the Optimal Control Problem 
%
% Double Integrator Tracking Problem
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

[problem,guess]=DoubleIntegratorTracking;          % Fetch the problem definition
options= problem.settings(30);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );

%% figure
xx=linspace(solution.T(1,1),solution.tf,1000);


figure
hold on
plot(xx,speval(solution,'X',1,xx),'r-' )
plot(xx,speval(solution,'X',2,xx),'b-' )
plot(tv,xv(:,1),'k-.' )
plot(tv,xv(:,2),'k-.' )
xlabel('Time [s]')
ylabel('States')
legend('Position [m]','Velocity [m/s]')
grid on

figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u [N]')

figure
plot(xx,(speval(solution,'X',1,xx)-5.*sin(xx')))
hold on
plot(tv,(xv(:,1)-5.*sin(tv)),'k-.' )
xlabel('Time [s]')
grid on
ylabel('Tracking error')
legend('Position error [m]')

