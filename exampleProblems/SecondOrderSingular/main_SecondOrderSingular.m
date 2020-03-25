% Main script to solve the Optimal Control Problem 
%
% Second order singular regulator problem
%
% The problem was originally presented by: 
% G.M. Aly. 1978.   The computation of optimal singular control.Internat. J. Control28, 5 (1978), 681–688. 
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
[problem,guess]=SecondOrderSingular;          % Fetch the problem definition
options= problem.settings(100);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);


%% figure
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
hold on
plot(xx,speval(solution,'X',1,xx),'b-' )
plot(xx,speval(solution,'U',1,xx),'r-' )
xlabel('Time [s]')
ylabel('States')
legend('x1 [-]','x2 [-]')
grid on

figure
hold on
if length(solution.T)==length(solution.dU)
    plot(solution.T(:,1),solution.dU(:,1),'b-' );
else
    plot(solution.T(1:end-1,1),solution.dU(:,1),'b-' );
end
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u [N]')
    