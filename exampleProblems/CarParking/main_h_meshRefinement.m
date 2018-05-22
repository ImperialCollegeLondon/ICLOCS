% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Minimum Time Parallel Parking
% The problem was adapted from 
% B. Li, K. Wang, and Z. Shao, "Time-optimal maneuver planning in automatic parallel parking using a simultaneous dynamic optimization approach". IEEE Transactions on Intelligent Transportation Systems, 17(11), pp.3263-3274, 2016.
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

options= settings_h(40);                  % Get options and solver settings 
[problem,guess]=CarParking;          % Fetch the problem definition
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
posx=solution.X(:,1);
posy=solution.X(:,2);
theta=solution.X(:,4);

A_x=posx+(data.data.auxdata.l_axes+data.data.auxdata.l_front).*cos(theta)-data.data.auxdata.b_width.*sin(theta);
B_x=posx+(data.data.auxdata.l_axes+data.data.auxdata.l_front).*cos(theta)+data.data.auxdata.b_width.*sin(theta);
C_x=posx-data.data.auxdata.l_rear.*cos(theta)+data.data.auxdata.b_width.*sin(theta);
D_x=posx-data.data.auxdata.l_rear.*cos(theta)-data.data.auxdata.b_width.*sin(theta);

A_y=posy+(data.data.auxdata.l_axes+data.data.auxdata.l_front).*sin(theta)+data.data.auxdata.b_width.*cos(theta);
B_y=posy+(data.data.auxdata.l_axes+data.data.auxdata.l_front).*sin(theta)-data.data.auxdata.b_width.*cos(theta);
C_y=posy-data.data.auxdata.l_rear.*sin(theta)-data.data.auxdata.b_width.*cos(theta);
D_y=posy-data.data.auxdata.l_rear.*sin(theta)+data.data.auxdata.b_width.*cos(theta);


%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
xv=linspace(problem.states.xl(1),problem.states.xu(1),10000);
hold on
plot([problem.states.xl(1),problem.states.xu(1)],[data.data.auxdata.CL data.data.auxdata.CL],'r-')
plot(xv,(-tanh(50*(xv-0.1))+tanh(50*(xv-data.data.auxdata.SL+0.1)))/2.*data.data.auxdata.SW,'r-')
for i=1:1:size(A_x,1)
    plot([A_x(i) B_x(i) C_x(i) D_x(i) A_x(i)], [A_y(i) B_y(i) C_y(i) D_y(i) A_y(i)], 'b-')
end
plot(solution.X(:,1),solution.X(:,2),'ko-','linewidth', 2)
xlim([-2 10])
xlabel('Position x [m]')
ylabel('Position y [m]')
grid on


figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,1,xx),'b-' )
xlabel('Position x [m]')
ylabel('Time [s]')
grid on


figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,2,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,2,xx),'b-' )
xlabel('Position y [m]')
ylabel('Time [s]')
grid on


%%
figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,3,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,3,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.states.xl(3), problem.states.xl(3)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(3), problem.states.xu(3)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
% 
figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T(:,1)),'bo' )
hold on
plot(xx,speval(solution.Up,1,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
grid on

% 
figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,4,xx)*180/pi,'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Orientation Angle [deg]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,5,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,5,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.states.xl(5), problem.states.xl(5)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(5), problem.states.xu(5)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Steering Angle [deg]')
grid on

figure
plot(solution.T(1:end-1,1),solution.dU(:,1)*180/pi,'-o')
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.url(1), problem.inputs.url(1)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uru(1), problem.inputs.uru(1)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Jerk [m/s^3]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,2,solution.T)*180/pi,'-o')
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Control Input (steering rate) [deg/s]')
grid on