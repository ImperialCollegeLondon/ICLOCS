% MAIN - Main script to solve the Optimal Control Problem
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all;close all;format compact;
global sol;  
sol=[];                             % Initialize solution structure

options= settings_hp(15,4);                  % Get options and solver settings 
[problem,guess]=F50MinTimeFlight;          % Fetch the problem definition

errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;

problem_iter=problem;
while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
    if i~=1
        [ problem_iter,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, i );
    end
    [infoNLP,data,options]=transcribeOCP(problem_iter,guess,options); % Format for NLP solver
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

%%
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
ylim([-1 9]*10^5)

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
text(8.52e05, 7.8e05,'NO FLIGHT ZONE','Color','white','FontSize',14)
xlim([8.5 9.1]*10^5)
ylim([7.5 8.1]*10^5)
%%
figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,1,solution.TSeg_Bar,xx),'b-' )
% plot([solution.T(1,1); solution.tf],[problem.states.xl(1), problem.states.xl(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(1), problem.states.xu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,4,solution.TSeg_Bar,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('True Airspeed [m/s]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,5,solution.TSeg_Bar,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,5,solution.TSeg_Bar,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on


figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,6,solution.TSeg_Bar,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,6,solution.TSeg_Bar,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Tracking Angle [deg]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,7,solution.TSeg_Bar,[solution.T(:,1); solution.tf])/9.81,'bo' )
hold on
plot(xx,speval(solution.Xp,7,solution.TSeg_Bar,xx)/9.81,'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,1,solution.TSeg_Bar,solution.T)*180/pi,'bo')
hold on
plot(xx,speval(solution.Up,1,solution.TSeg_Bar,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Angle of attack [deg]')
grid on


figure
plot(solution.T(:,1),speval(solution.Up,2,solution.TSeg_Bar,solution.T)*180/pi,'bo')
hold on
plot(xx,speval(solution.Up,2,solution.TSeg_Bar,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Roll angle [deg]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,3,solution.TSeg_Bar,solution.T),'bo')
hold on
plot(xx,speval(solution.Up,3,solution.TSeg_Bar,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(3), problem.inputs.ul(3)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(3), problem.inputs.uu(3)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Throttle Setting [-]')
grid on

