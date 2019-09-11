% Main script to solve the Optimal Control Problem
%
% OrbitRaising - Orbit Raising Problem
% 
% The problem was originally presented by: 
% A. E. Bryson, and Y.C. Ho, "Applied Optimal Control", Hemisphere Publishing, New York, 1975.
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

[problem,guess]=OrbitRaising;          % Fetch the problem definition
% options= problem.settings(4,4);              % for hp method
options= problem.settings(40);                  % for h method
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.01 );


%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
plot(solution.T,speval(solution,'X',1,solution.T),'ro' )
hold on
plot(solution.T,speval(solution,'X',2,solution.T),'bo' )
plot(solution.T,speval(solution,'X',3,solution.T),'mo' )
plot(solution.T,speval(solution,'X',4,solution.T),'go' )
plot(xx,speval(solution,'X',1,xx),'r-' )
plot(xx,speval(solution,'X',2,xx),'b-' )
plot(xx,speval(solution,'X',3,xx),'m-' )
plot(xx,speval(solution,'X',4,xx),'g-' )
plot(tv,xv(:,1),'k-.')
plot(tv,xv(:,2),'k-.')
plot(tv,xv(:,3),'k-.')
plot(tv,xv(:,4),'k-.')
xlabel('Time')
ylim([0 2.5])
ylabel('States')
legend('r','theta','V_r','V_{theta}')
grid on


figure
plot(solution.T(:,1),speval(solution,'U',1,solution.T),'ro' )
hold on
plot(solution.T(:,1),speval(solution,'U',2,solution.T),'bo' )
plot(xx,speval(solution,'U',1,xx),'r-' )
plot(xx,speval(solution,'U',2,xx),'b-' )
plot(tv,uv(:,1),'k-.')
plot(tv,uv(:,2),'k-.')
xlabel('Time')
grid on
ylabel('Control Inputs')
legend('u1','u2')




