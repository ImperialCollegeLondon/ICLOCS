% main_Brachistochrone - Main script to solve the Optimal Control Problem
%
% Brachistochrone Problem
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




%% Solve the problem
clear all;close all;format compact;
[problem,guess]=Brachistochrone;          % Fetch the problem definition
options= problem.settings(100);           % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);

%% Generate Figures
xx=linspace(solution.T(1,1),solution.tf,100000);

figure(100)
hold on
plot(speval(solution,'X',1,xx),speval(solution,'X',2,xx))
xlabel('X [m]')
ylabel('Y [m]')
grid on
set(gca,'YDir','reverse');

figure(101)
hold on
plot(xx,speval(solution,'X',3,xx))
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on

figure(102)
hold on
plot(xx,speval(solution,'U',1,xx))
xlabel('Time [s]')
ylabel('Control [rad]')
grid on