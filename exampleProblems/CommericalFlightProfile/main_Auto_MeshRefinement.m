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

options= settings_Auto(40);                  % Get options and solver settings 
[problem,guess]=F50MinTimeFlight;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
iterHistory=zeros(1,2);

solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;

problem_iter=problem;
while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax   
%     if i~=1
%         [ problem_iter,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, i );
%     end
    [infoNLP,data,options]=transcribeOCP(problem_iter,guess,options); % Format for NLP solver
    [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
    [solution]=output(problem,solution,options,data,4);         % Output solutions
    
    
    maxAbsError=max(abs(solution.Error));
    maxAbsConstraintError=max(solution.ConstraintError);
    errorHistory(i,:)=maxAbsError;
    iterHistory(i)=status.iter;
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
        xx=linspace(solution.T(1,1),solution.tf,1000);

        figure
        center=[data.data.auxdata.obs_epos_1 data.data.auxdata.obs_npos_1];
        obspos = [center-data.data.auxdata.obs_r_1 2*data.data.auxdata.obs_r_1 2*data.data.auxdata.obs_r_1];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
        hold on
        center=[data.data.auxdata.obs_epos_2 data.data.auxdata.obs_npos_2];
        obspos = [center-data.data.auxdata.obs_r_2 2*data.data.auxdata.obs_r_2 2*data.data.auxdata.obs_r_2];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
        center=[data.data.auxdata.obs_epos_3 data.data.auxdata.obs_npos_3];
        obspos = [center-data.data.auxdata.obs_r_3 2*data.data.auxdata.obs_r_3 2*data.data.auxdata.obs_r_3];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
        center=[data.data.auxdata.obs_epos_4 data.data.auxdata.obs_npos_4];
        obspos = [center-data.data.auxdata.obs_r_4 2*data.data.auxdata.obs_r_4 2*data.data.auxdata.obs_r_4];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
        center=[data.data.auxdata.obs_epos_5 data.data.auxdata.obs_npos_5];
        obspos = [center-data.data.auxdata.obs_r_5 2*data.data.auxdata.obs_r_5 2*data.data.auxdata.obs_r_5];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');

        plot(solution.X(:,3),solution.X(:,2),'b-','LineWidth',2)
        xlabel('East Position [m]')
        ylabel('North Position [m]')
        grid on
        plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
        plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
        text(problem.states.x0(3)+20000,problem.states.x0(2),'ORG')
        text(problem.states.xfl(3)+10000,problem.states.xfl(2)+10000,'DES')
        xlim([-1 10]*10^5)
        ylim([-0.5 9]*10^5)

        figure
        center=[data.data.auxdata.obs_epos_1 data.data.auxdata.obs_npos_1];
        obspos = [center-data.data.auxdata.obs_r_1 2*data.data.auxdata.obs_r_1 2*data.data.auxdata.obs_r_1];
        rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
        hold on
        plot(solution.X(:,3),solution.X(:,2),'b-','LineWidth',2)
        xlabel('East Position [m]')
        ylabel('North Position [m]')
        grid on
        plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
        plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
        text(problem.states.x0(3)+10000,problem.states.x0(2),'ORG')
        text(problem.states.xfl(3)-10000,problem.states.xfl(2)-5000,'DES')
        % text(8.25e05, 7.45e05,'NO FLIGHT ZONE','Color','white','FontSize',14)
        xlim([8.5 9.1]*10^5)
        ylim([7.5 8.1]*10^5)

        %%
        figure
        subplot(2,1,1)
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,1,xx),'b-' )
        % plot([solution.T(1,1); solution.tf],[problem.states.xl(1), problem.states.xl(1)],'r-' )
        plot([solution.T(1,1); solution.tf],[problem.states.xu(1), problem.states.xu(1)],'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Altitude [m]')
        grid on

        subplot(2,1,2)
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf]),'bo' )
        hold on
        plot(xx,speval(solution.Xp,4,xx),'b-' )
        plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)],'r-' )
        plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)],'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('True Airspeed [m/s]')
        grid on

        figure
        subplot(2,1,1)
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,5,[solution.T(:,1); solution.tf])*180/pi,'bo' )
        hold on
        plot(xx,speval(solution.Xp,5,xx)*180/pi,'b-' )
        hold on
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Flight Path Angle [deg]')
        grid on


        subplot(2,1,2)
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,6,[solution.T(:,1); solution.tf])*180/pi,'bo' )
        hold on
        plot(xx,speval(solution.Xp,6,xx)*180/pi,'b-' )
        hold on
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Tracking Angle [deg]')
        grid on

        figure
        subplot(2,1,1)
        plot([solution.T(:,1); solution.tf],speval(solution.Xp,7,[solution.T(:,1); solution.tf])/9.81,'bo' )
        hold on
        plot(xx,speval(solution.Xp,7,xx)/9.81,'b-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Aircraft Mass [kg]')
        grid on

        subplot(2,1,2)
        plot(solution.T(:,1),speval(solution.Up,1,solution.T)*180/pi,'bo')
        hold on
        plot(xx,speval(solution.Up,1,xx)*180/pi,'b-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)]*180/pi,'r-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)]*180/pi,'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Angle of attack (Control) [deg]')
        grid on


        figure
        subplot(2,1,1)
        plot(solution.T(:,1),speval(solution.Up,2,solution.T)*180/pi,'bo')
        hold on
        plot(xx,speval(solution.Up,2,xx)*180/pi,'b-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Roll angle (Control) [deg]')
        grid on

        subplot(2,1,2)
        % plot(solution.T(:,1),sqrt(speval(solution.Up,3,solution.T)),'bo')
        % hold on
        % plot(xx,sqrt(speval(solution.Up,3,xx)),'b-' )
        plot(solution.T(:,1),speval(solution.Up,3,solution.T),'bo')
        hold on
        plot(xx,speval(solution.Up,3,xx),'b-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.ul(3), problem.inputs.ul(3)],'r-' )
        plot([solution.T(1,1); solution.tf],[problem.inputs.uu(3), problem.inputs.uu(3)],'r-' )
        xlim([0 solution.tf])
        xlabel('Time [s]')
        ylabel('Throttle Setting (Control) [-]')
        grid on
end

