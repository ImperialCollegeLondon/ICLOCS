% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% BangBang Control (Double Integrator Minimum Time Repositioning) Problem
%
% The problem was adapted from Example 4.11 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all;close all;format compact;


[problem,guess]=ShuttleReentryTrajectory;          % Fetch the problem definition
options= problem.settings(6,4);
% options= problem.settings(30);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.1 );


%% figure
xx=linspace(solution.T(1,1),solution.tf,100000);

figure
hold on
plot(xx,(speval(solution,'X',1,xx))/100000,'b-' )
xlim([0 solution.tf])
plot(tv,(xv(:,1))/100000,'k-.')
xlabel('Time [s]')
ylabel('h [100,000 ft]')
grid on


figure
hold on
plot(xx,speval(solution,'X',2,xx)*180/pi,'b-' )
xlim([0 solution.tf])
plot(tv,xv(:,2)*180/pi,'k-.')
xlabel('Time [s]')
ylabel('Longitude [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',3,xx)*180/pi,'b-' )
plot(tv,xv(:,3)*180/pi,'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Latitude [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',4,xx)/1000,'b-' )
xlim([0 solution.tf])
plot(tv,xv(:,4)/1000,'k-.')
xlabel('Time [s]')
ylabel('Velocity [1000 ft/s]')
grid on

figure
hold on
plot(xx,speval(solution,'X',5,xx)*180/pi,'b-' )
plot(tv,xv(:,5)*180/pi,'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',6,xx)*180/pi,'b-' )
plot(tv,xv(:,6)*180/pi,'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Azimuth [deg]')
grid on

%%
figure
hold on
plot(xx,speval(solution,'U',1,xx)*180/pi,'b-' )
xlim([0 solution.tf])
plot(tv,uv(:,1)*180/pi,'k-.')
xlabel('Time [s]')
ylabel('Angle of Attack [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'U',2,xx)*180/pi,'b-' )
plot(tv,uv(:,2)*180/pi,'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Bank Angle [deg]')
grid on
