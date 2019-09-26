% main_GoddardRocketThreePhase - Main script to solve the Optimal Control Problem with a multi-phase formulation
%
% Goddard Rocket Problem (Three-phase Formulation)
%
% The problem was adapted from Example 4.9 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
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

options.mp= settings_GoddardRocketThreePhase;                  % Get options and solver settings 
[problem,guess,options.phaseoptions]=GoddardRocketThreePhase;          % Fetch the problem definition
[solution,MRHistory]=solveMyProblem( problem,guess,options);

%%
for i=1:length(solution.phaseSol)
    sol=solution.phaseSol{i};
    xx=linspace(sol.t0,sol.tf,1000);
    
    
    figure(100)
    hold on
    plot(xx,speval(sol,'X',1,xx),'linewidth',2)
    xlabel('Time [s]')
    ylabel('Altitude [ft]')
    grid on
    
    figure(101)
    hold on
    plot(xx,speval(sol,'X',2,xx),'linewidth',2)
    xlabel('Time [s]')
    ylabel('Velocity [ft/s]')
    grid on
    
    figure(102)
    hold on
    plot(xx,speval(sol,'X',3,xx),'linewidth',2)
    xlabel('Time [s]')
    ylabel('Mass [lbm]')
    grid on
    
    figure(103)
    hold on
    plot(xx,speval(sol,'U',1,xx),'linewidth',2)
    ylim([problem.phases{3}.inputs.ul problem.phases{1}.inputs.uu])
    xlabel('Time [s]')
    grid on
    ylabel('Control Input (Thrust) [lbf]')

end



    