function [problem,guess,phaseoptions] = myMultiPhaseProblem
%myMultiPhaseProblem - Template file for optimal control problem definition using multi-phase formulation
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%    phaseoptions - options for each phases
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


% Initial and final time for different phases. Let t_min(end)=t_max(end) if tf is fixed.
problem.mp.time.t_min=[t0_min t1_min ...];     
problem.mp.time.t_max=[t0_max t1_max ...];  
guess.mp.time=[t0_guess t1_guess ...];

% Parameters bounds. pl=< p <=pu
problem.mp.parameters.pl=[];
problem.mp.parameters.pu=[];
guess.mp.parameters=[];

% Bounds for linkage boundary constraints bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
problem.mp.constraints.bll.linear=[bl1_linear_lowerbound bl2_linear_lowerbound ...];
problem.mp.constraints.blu.linear=[bl1_linear_upperbound bl2_linear_upperbound ...];
problem.mp.constraints.blTol.linear=[eps_bl1_linear_bounds eps_bl2_linear_bounds ...]; 

problem.mp.constraints.bll.nonlinear=[bl1_nonlinear_lowerbound bl2_nonlinear_lowerbound ...];
problem.mp.constraints.blu.nonlinear=[bl1_nonlinear_upperbound bl2_nonlinear_upperbound ...];
problem.mp.constraints.blTol.nonlinear=[eps_bl1_nonlinear_bounds eps_bl2_nonlinear_bounds ...]; 

% Get function handles
problem.mp.linkfunctions=@bclink;

% Store the necessary problem parameters used in the functions
problem.mp.data = [];

% Define different phases of OCP
[problem.phases{1},guess.phases{1}] = myMultiPhaseProblem_Phase1(problem.mp, guess.mp);
[problem.phases{2},guess.phases{2}] = myMultiPhaseProblem_Phase2(problem.mp, guess.mp);
...

% Each phase could use different discretization method
phaseoptions{1}=problem.phases{1}.settings(Nps); % h method for phase 1
phaseoptions{2}=problem.phases{2}.settings(Nps,Npd); % hp method for phase 2
...
%------------- END OF CODE --------------


function [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)

% bclink - Returns the evaluation of the linkage boundary constraints: bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
%
% Syntax:  [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    vdat- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    blc_linear - column vector containing the evaluation of the linear linkage boundary constraint functions
%    blc_nonlinear - column vector containing the evaluation of the nonlinear linkage boundary constraint functions
%
%------------- BEGIN CODE --------------

% Variable of different phase could be called using syntex of, for example, xf{n_phase}(n): the value of nth state at tf of phase number 'n_phase'

% linear linkage constraints 
blc_linear(1,:)=blc_linear_1(x0,xf,u0,uf,p,t0,tf,vdat);
blc_linear(2,:)=blc_linear_2(x0,xf,u0,uf,p,t0,tf,vdat);
...

% nonlinear linkage constraints 
blc_nonlinear(1,:)=blc_nonlinear_1(x0,xf,u0,uf,p,t0,tf,vdat);
blc_nonlinear(2,:)=blc_nonlinear_2(x0,xf,u0,uf,p,t0,tf,vdat);
...
%------------- END OF CODE --------------
