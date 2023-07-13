% Main script to solve the Optimal Control Problem 
%
% Hypersensitive problem
%
% The problem was originally presented by: 
% A.V. Rao, and K.D. Mease, "Eigenvector approximate dichotomic basis method for solving hyper?sensitive optimal control problems", Optimal Control Applications and Methods, 21(1), pp.1-19, 2000
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

[problem,guess]=Hypersensitive;          % Fetch the problem definition
options= problem.settings(1,8);                  % for hp method
% options= problem.settings(40);                  % for h method
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );

%% figure
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
plot(xx,speval(solution,'X',1,xx),'b-' )
hold on
plot(tv,xv(:,1),'k-.' )
xlabel('Time [s]')
ylabel('x(t)')
grid on

figure
hold on
plot(xx,speval(solution,'U',1,xx),'b-' )
plot(tv,uv(:,1),'k-.' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
xlabel('Time [s]')
grid on
ylabel('u(t)')
    

