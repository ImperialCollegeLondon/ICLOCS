% Main script to solve the Optimal Control Problem 
%
% Aircraft go around in the present of wind-shear problem
%
% The problem was originally presented by: 
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 1: Necessary conditions", Journal of Optimization Theory and Applications, 70(1), pp 1-23, 1991.
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 2: Multiple shooting and homotopy", Journal of Optimization Theory and Applications, 70(2), pp 223-254, 1991.
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% This implementation contains modifications
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

[problem,guess]=WindshearGoAround;          % Fetch the problem definition
% options= problem.settings(5,8);                  % for hp method
options= problem.settings(40);                  % for h method
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% Generation of Plots
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
plot(solution.X(:,1),solution.X(:,2),'bo-')
hold on
h1=plot([solution.X(1,1); solution.X(end,1)],[solution.p, solution.p],'k--' );
xlabel('Position [ft]');
ylabel('Altitude [ft]');
xlim([solution.X(1,1); solution.X(end,1)])
legend(h1,'minimum altitude')
grid on

figure
hold on
plot(xx,speval(solution,'X',3,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity [ft/s]')
grid on


figure
hold on
plot(xx,speval(solution,'X',4,xx)*180/pi,'b-' )
hold on
plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
hold on
plot(xx,speval(solution,'X',1,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Position [ft]')
grid on


figure
hold on
plot(xx,speval(solution,'U',1,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu]*180/pi,'r--' )
xlim([0 solution.tf])
xlabel('Time [s]');
ylabel('Angle of attack [deg]');
grid on

figure
if length(solution.T)==length(solution.dU)
    plot(solution.T(1:end,1),solution.dU*180/pi,'-b')
else
    plot(solution.T(1:end-1,1),solution.dU*180/pi,'-b')
end
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.url, problem.inputs.url]*180/pi,'r--' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uru, problem.inputs.uru]*180/pi,'r--' )
xlim([0 solution.tf])
xlabel('Time [s]');
ylabel('Angle of attack rate [deg/s]');
grid on

