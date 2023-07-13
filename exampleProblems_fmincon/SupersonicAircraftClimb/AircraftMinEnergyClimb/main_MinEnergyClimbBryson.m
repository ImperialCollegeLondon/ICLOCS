% Main script to solve the Optimal Control Problem
%
% Supersonic Aircraft Minimum Time-to-climb 
%
% The problem was adapted from the supersonic aircraft minimum time-to-climb problem originally presented by
% A. E. Bryson, M. N. Desai, and W. C. Hoffman, "Energy-State Approximation in Performance Optimization of Supersonic Aircraft," Journal of Aircraft, Vol. 6, No. 6, November-December, 1969, pp. 481-488. 
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% With aerodynamic data and modifications to thrust data by:
% M.A. Patterson and A.V. Rao, User's Manual, "GPOPS-II: A General Purpose MATLAB Software for Solving Multiple-Phase Optimal Control Problems, 2.3 edition", 2016
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

[problem,guess]=MinEnergyClimbBryson;          % Fetch the problem definition
% options= problem.settings(10,4);              % for hp method
options= problem.settings(20);                  % for h method
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figure
xx=linspace(solution.T(1,1),solution.tf,10000);

figure
plot(speval(solution,'X',2,xx)/100,speval(solution,'X',1,xx),'b-')
hold on
plot(xv(:,2)/100,xv(:,1),'k-.')
ylim([0 20000])
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
grid on

figure
hold on
plot(xx,speval(solution,'X',3,xx)*180/pi,'b-' )
plot(tv,xv(:,3)*180/pi,'k-.')
plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',1,xx),'b-' )
plot(tv,xv(:,1),'k-.')
xlim([0 solution.tf])
ylim([0 20000])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure
hold on
plot(xx,speval(solution,'X',2,xx)/100,'b-' )
plot(tv,xv(:,2)/100,'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
grid on

figure
hold on
plot(xx,speval(solution,'X',4,xx),'b-' )
plot(tv,xv(:,4),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
grid on

