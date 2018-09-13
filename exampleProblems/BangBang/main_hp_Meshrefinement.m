% main_hp_MeshRefinement - Main script to solve the Optimal Control Problem with hp-typed mesh and refinement
%
% BangBang Control (Double Integrator Minimum Time Repositioning) Problem
%
% The problem was adapted from Example 4.11 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
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

options= settings_hp(1,5);                  % Get options and solver settings 
[problem,guess]=BangBang;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
iterHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;
minItervalScale=1;
while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
   
    [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
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
MeshRefinementHistory.iterHistory=iterHistory;
MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;

%%
xx=linspace(solution.T(1,1),solution.tf,10000);

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'ro' )
hold on
plot([solution.T(:,1); solution.tf],speval(solution.Xp,2,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bd' )
plot(xx,speval(solution.Xp,1,solution.TSeg_Bar,xx),'r-' )
plot(xx,speval(solution.Xp,2,solution.TSeg_Bar,xx),'b-.' )
xlabel('Time [s]')
ylabel('States')
legend('Position [m]','Velocity [m/s]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Up,1,solution.TSeg_Bar,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul, problem.inputs.ul],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu],'r-' )
plot(xx,speval(solution.Up,1,solution.TSeg_Bar,xx),'b-' )
xlabel('Time [s]')
grid on
ylabel('Control Input')
legend('u [N]')
 