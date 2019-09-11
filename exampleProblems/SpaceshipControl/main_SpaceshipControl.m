% Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Optimal Control for a Spaceship
%
% The problem was originally presented by: 
% Battin, R.H.: An Introduction to the Mathematics and Methods of Astrodynamics, Revised Edition. AIAA Education Series (1999)
% This implementation was adapted from
% Knauer M., Büskens C. (2019) Real-Time Optimal Control Using TransWORHP and WORHP Zen. In: Fasano G., Pintér J. (eds) Modeling and Optimization in Space Engineering. Springer Optimization and Its Applications, vol 144. Springer, Cham
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

options= settings_SpaceshipControl(6,3);
% options= settings_SpaceshipControl(20);                  % Get options and solver settings 
[problem,guess]=SpaceshipControl;          % Fetch the problem definition
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figure
xx=linspace(solution.T(1,1),solution.tf,100000);

figure
plot(xx,speval(solution,'X',1,xx),'b-' )
hold on
plot(tv,xv(:,1),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Position x')
grid on

figure
plot(xx,speval(solution,'X',2,xx),'b-' )
hold on
plot(tv,xv(:,2),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Position y')
grid on

figure
plot(xx,speval(solution,'X',3,xx),'b-' )
hold on
plot(tv,xv(:,3),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity x')
grid on

figure
plot(xx,speval(solution,'X',4,xx),'b-' )
hold on
plot(tv,xv(:,4),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity y')
grid on


% 
figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control x')
grid on

figure
hold on
plot(xx,speval(solution,'U',2,xx),'b-' )
plot(tv,uv(:,2),'k-.')
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control y')
grid on

%%