% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Fed-batch fermentor optimal control problem
%
% The problem was adapted from Example 2, Section 12.4.2 of
% J.E. Cuthrell and L. T. Biegler. Simultaneous optimization and solution methods for batch reactor control profiles. Computers and Chemical Engineering, 13:49-62, 1989.
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
[problem,guess]=BatchFermentor;          % Fetch the problem definition
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
xx=linspace(solution.T(1,1),solution.T(end,1),1000);
figure
plot(solution.T(:,1),speval(solution.Xp,1,solution.T(:,1)),'ro' )
hold on
plot(solution.T(:,1),speval(solution.Xp,2,solution.T(:,1)),'bo' )
plot(solution.T(:,1),speval(solution.Xp,3,solution.T(:,1)),'mo' )
plot(solution.T(:,1),speval(solution.Xp,4,solution.T(:,1)),'ko' )
plot(xx,speval(solution.Xp,1,xx),'r-' )
plot(xx,speval(solution.Xp,2,xx),'b-' )
plot(xx,speval(solution.Xp,3,xx),'m-' )
plot(xx,speval(solution.Xp,4,xx),'k-' )
xlim([0 solution.tf])
xlabel('Time [hrs]')
ylabel('States')
legend('x [g/l]','p [g/l]','s [g/l]','v [l]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-bo' )
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [hrs]')
grid on
ylabel('Control Input')
legend('u [g/hrs]')
   
    