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
% global A1 B1 w1 W1 tau1 M1
global sol;  
sol=[];                             % Initialize solution structure

options= settings_hp(1,8);                  % Get options and solver settings 
[problem,guess]=F50MinTimeFlight;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;
minItervalScale=1;
while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
    
    [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
    [solution,infoNLP1,data] = solveNLP(infoNLP,data);      % Solve the NLP
    [solution]=output(problem,solution,options,data,4);         % Output solutions

    
    maxAbsError=max(abs(solution.Error));
    maxAbsConstraintError=max(solution.MaxConstraintError);
    errorHistory(i,:)=maxAbsError;
    ConstraintErrorHistory(i,:)=maxAbsConstraintError;
    timeHistory(i)=solution.computation_time;
    solutionHistory{i}=solution;
    i=i+1;
    %%
    if i~=1 && all(maxAbsError==errorHistory(i-1,:))
        minItervalScale=minItervalScale*0.9;
    end
    
    if (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
        [ options, guess ] = doMeshRefinement( options, problem, guess, data, solution, minItervalScale );
    end

end

MeshRefinementHistory.errorHistory=errorHistory;
MeshRefinementHistory.timeHistory=timeHistory;
MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;

%%
sos = ppval(problem.data.auxdata.us1976sos,solution.X(:,1));
Mach = solution.X(:,2)./sos;

%%
xx=linspace(solution.T(1,1),solution.tf,1000);
t0=solution.z(end-data.sizes{1}+1);tf=solution.z(end);
if data.options.adaptseg==1 
    TSeg_Bar=solution.z((end-nt+1):end);
else
    TSeg_Bar=(tf-t0)/2*data.tau_segment'+(tf+t0)/2; %Time at start/end of each segmen
end

figure
plot(Mach,solution.X(:,1),'bo-')
xlabel('Airspeed [Mach Number]')
ylabel('Altitude [m]')
grid on

figure
plot(ppvalL(solution.Xp,2,TSeg_Bar,[solution.T(:,1); solution.tf])/100,ppvalL(solution.Xp,1,TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(ppvalL(solution.Xp,2,TSeg_Bar,xx)/100,ppvalL(solution.Xp,1,TSeg_Bar,xx),'b-' )
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
grid on

figure
plot([solution.T(:,1); solution.tf],ppvalL(solution.Xp,3,TSeg_Bar,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,ppvalL(solution.Xp,3,TSeg_Bar,xx)*180/pi,'b-' )
% st=length(solution.z)-data.sizes{1}+1;
% ed=length(solution.z);
% for i=st:ed
%     plot([solution.z(i) solution.z(i)],[min(solution.X(:,3)) max(solution.X(:,3))]*180/pi,'k--')
% end
hold on
plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
plot([solution.T(:,1); solution.tf],ppvalL(solution.Xp,1,TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,ppvalL(solution.Xp,1,TSeg_Bar,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure
plot([solution.T(:,1); solution.tf],ppvalL(solution.Xp,4,TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,ppvalL(solution.Xp,4,TSeg_Bar,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure
plot([solution.T(:,1); solution.tf],ppvalL(solution.Xp,2,TSeg_Bar,[solution.T(:,1); solution.tf])/100,'bo' )
hold on
plot(xx,ppvalL(solution.Xp,2,TSeg_Bar,xx)/100,'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
grid on

figure
plot([solution.T(:,1); solution.tf],ppvalL(solution.Up,1,TSeg_Bar,[solution.T(:,1); solution.tf])/100,'bo' )
hold on
plot(xx,ppvalL(solution.Up,1,TSeg_Bar,xx)/100,'b-' )
hold on
% plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul]*180/pi,'r-' )
% plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
grid on
%%

figure
plot([solution.T(:,1); solution.tf],solution.multipliers.lambda(:,1),'bo-')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Costate for Altitude')
grid on


figure
plot([solution.T(:,1); solution.tf],solution.multipliers.lambda(:,2),'bo-')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Costate for Airspeed')
grid on

figure
plot([solution.T(:,1); solution.tf],solution.multipliers.lambda(:,3),'bo-')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Costate for Flight Path Angle')
grid on

figure
plot([solution.T(:,1); solution.tf],solution.multipliers.lambda(:,4),'bo-')
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Costate for Mass')
grid on


