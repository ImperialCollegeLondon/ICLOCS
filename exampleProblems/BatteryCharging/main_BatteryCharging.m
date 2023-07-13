% Main script to solve the Optimal Control Problem 
%
% Battery Charging Problem
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

[problem,guess]=BatteryCharging;          % Fetch the problem definition
options= problem.settings(50);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113');

%% figure
tt=solution.T;
x1=speval(solution,'X',1,tt);
x2=speval(solution,'X',2,tt);
u1=speval(solution,'U',1,tt);
y=3.64+0.55*x1-0.72*x1.^2+0.75*x1.^3+x2+problem.data.R0*u1;

figure
subplot(2,2,1)
hold on
plot(tt,x1,'b-' ,'LineWidth',2)
plot([solution.T(1,1); solution.tf],[problem.states.xl(1), problem.states.xl(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(1), problem.states.xu(1)],'r-' )
xlabel('Time [s]')
ylabel('States: State-of-Charge')
grid on

% figure
subplot(2,2,2)
hold on
plot(tt,x2,'b-' ,'LineWidth',2)
plot([solution.T(1,1); solution.tf],[problem.states.xl(2), problem.states.xl(2)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(2), problem.states.xu(2)],'r-' )
xlabel('Time [s]')
ylabel('States: RC Voltage [V]')
grid on

% figure
subplot(2,2,3)
hold on
plot(tt,u1,'b-' ,'LineWidth',2)
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input: Current [A]')

% figure
subplot(2,2,4)
hold on
plot(tt,y,'b-' ,'LineWidth',2)
plot([solution.T(1,1); solution.tf],[problem.constraints.gl(1), problem.constraints.gl(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.constraints.gu(1), problem.constraints.gu(1)],'r-' )
xlabel('Time [s]')
ylabel('Output: voltage [V]')
grid on
