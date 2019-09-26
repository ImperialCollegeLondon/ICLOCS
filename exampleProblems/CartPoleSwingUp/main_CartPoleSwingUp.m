% main_h_MeshRefinement - Main script to solve the Cart Pole Swing-up Problem
%
% Cart Pole Swing-up optimal control problem
%
% The problem was adapted from Section 6 of
% Kelly M. An introduction to trajectory optimization: How to do your own direct collocation. SIAM Review. 2017 Nov 6;59(4):849-904.
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

[problem,guess]=CartPoleSwingUp;          % Fetch the problem definition
options= problem.settings(30);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);

problem.sim.functions=@CartPoleSwingUp_Dynamics_Sim_Exact;
[ tv1, xv1, uv1 ] = simulateSolution( problem, solution, 'ode23', 0.01 );

problem.sim.functions=@CartPoleSwingUp_Dynamics_Sim_Inexact;
[ tv2, xv2, uv2 ] = simulateSolution( problem, solution, 'ode23', 0.01 );




%% figures
xx=linspace(solution.T(1,1),solution.tf,100000);
figure
hold on
plot(xx,speval(solution,'X',3,xx),'linewidth',2)
plot(tv1,xv1(:,3),'k--','linewidth',2)
plot(tv2,xv2(:,3),'r--','linewidth',2)
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Position [m]')
legend('Optimization output','Simulation with exact dynamics','Simulation with modified dynamics')
grid on

figure
hold on
plot(xx,speval(solution,'X',4,xx),'linewidth',2)
plot(tv1,xv1(:,4),'k--','linewidth',2)
plot(tv2,xv2(:,4),'r--','linewidth',2)
xlim([0 solution.tf])
legend('Optimization output','Simulation with exact dynamics','Simulation with modified dynamics')
xlabel('Time [s]')
ylabel('Angle [rad]')
grid on


figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' ,'linewidth',2)
plot(tv1,uv1(:,1),'k--','linewidth',2)
plot(tv2,uv2(:,1),'r--','linewidth',2)
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
legend('Optimization output','Simulation with exact dynamics','Simulation with modified dynamics')
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input Force [N]')
