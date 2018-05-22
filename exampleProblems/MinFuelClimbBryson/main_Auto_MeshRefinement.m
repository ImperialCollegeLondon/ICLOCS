% main_Auto_MeshRefinement - Main script to solve the Optimal Control Problem with automatic mesh selection and refinement
%
% Supersonic Aircraft Minimum Fuel Climb 
% The problem was adapted from the supersonic aircraft minimum time-to-climb problem originally presented by 
% A. E. Bryson, M. N. Desai, and W. C. Hoffman, "Energy-State Approximation in Performance Optimization of Supersonic Aircraft," Journal of Aircraft, Vol. 6, No. 6, November-December, 1969, pp. 481-488. 
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% With aerodynamic data and modifications to thrust data by:
% M.A. Patterson and A.V. Rao, User's Manual, "GPOPS-II: A General Purpose MATLAB Software for Solving Multiple-Phase Optimal Control Problems, 2.3 edition", 2016
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%--------------------------------------------------------

clear all;close all;format compact;
global sol;  
sol=[];                             % Initialize solution structure

options= settings_Auto(20);                  % Get options and solver settings 
[problem,guess]=MinFuelClimbBryson;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;

while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax    
    [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
    [solution,infoNLP1,data] = solveNLP(infoNLP,data);      % Solve the NLP
    [solution]=output(problem,solution,options,data,4);         % Output solutions
    
    
    maxAbsError=max(abs(solution.Error));
    maxAbsConstraintError=max(solution.ConstraintError);
    errorHistory(i,:)=maxAbsError;
    ConstraintErrorHistory(i,:)=maxAbsConstraintError;
    timeHistory(i)=solution.computation_time;
    solutionHistory{i}=solution;
    
    if (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
        [ options, guess ] = doMeshRefinement( options, problem, guess, data, solution, i );
    end
    i=i+1;

end

MeshRefinementHistory.errorHistory=errorHistory;
MeshRefinementHistory.timeHistory=timeHistory;
MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;

%%
if (strcmp(options.transcription,'globalLGR')) || (strcmp(options.transcription,'hpLGR'))
        sos = ppval(problem.data.Atomssos,solution.X(:,1));
        Mach = solution.X(:,2)./sos;
        xx=linspace(solution.T(1,1),solution.tf,1000);

        figure
        plot(Mach,solution.X(:,1),'bo-')
        xlabel('Airspeed [Mach Number]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot(speval(solution.Xp,2,solution.TSeg_Bar,[solution.T(:,1); solution.tf])/100,speval(solution.Xp,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(speval(solution.Xp,2,solution.TSeg_Bar,xx)/100,speval(solution.Xp,1,solution.TSeg_Bar,xx),'b-' )
        xlabel('Airspeed [100 m/s]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,3,solution.TSeg_Bar,[solution.T(:,1); solution.tf])*180/pi,'bo' )
        hold on
        plot(xx,speval(solution.Xp,3,solution.TSeg_Bar,xx)*180/pi,'b-' )
        hold on
        plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
        plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Flight Path Angle [deg]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,1,solution.TSeg_Bar,xx),'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,4,solution.TSeg_Bar,xx),'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Aircraft Mass [kg]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,2,solution.TSeg_Bar,[solution.T(:,1); solution.tf])/100,'bo' )
        hold on
        plot(xx,speval(solution.Xp,2,solution.TSeg_Bar,xx)/100,'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Velocity [100 m/s]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Up,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf])/100,'bo' )
        hold on
        plot(xx,speval(solution.Up,1,solution.TSeg_Bar,xx)/100,'b-' )
        hold on
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Control Input (angle of attack) [deg]')
        grid on
else
        sos = ppval(problem.data.Atomssos,solution.X(:,1));
        Mach = solution.X(:,2)./sos;
        xx=linspace(solution.T(1,1),solution.tf,1000);

        figure
        plot(Mach,solution.X(:,1),'bo-')
        xlabel('Airspeed [Mach Number]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot(solution.X(:,2)/100,solution.X(:,1),'bo-')
        xlabel('Airspeed [100 m/s]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,3,[solution.T(:,1); solution.tf])*180/pi,'bo' )
        hold on
        plot(xx,speval(solution.Xp,3,xx)*180/pi,'b-' )
        hold on
        plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
        plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Flight Path Angle [deg]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,1,xx),'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Altitude [m]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,2,[solution.T(:,1); solution.tf])/100,'bo' )
        hold on
        plot(xx,speval(solution.Xp,2,xx)/100,'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Velocity [100 m/s]')
        grid on

        figure
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,4,xx),'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Aircraft Mass [kg]')
        grid on

        figure
        plot(solution.T(:,1),speval(solution.Up,1,solution.T)*180/pi,'-o')
        hold on
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Control Input (angle of attack) [deg]')
        grid on
end