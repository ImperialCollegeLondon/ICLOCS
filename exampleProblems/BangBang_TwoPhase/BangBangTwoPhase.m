function [problem,guess,phaseoptions] = BangBangTwoPhase
%BangBang - BangBang Control (Double Integrator Minimum Time Repositioning) Problem with a multi-phase formulation
%
% The problem was adapted from Example 4.11 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% MAT-files required: none
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
problem.mp.time.t_min=[0 0 0];     
problem.mp.time.t_max=[0 70 70]; 
guess.mp.time=[0 30 60];

% Parameters bounds. pl=< p <=pu
problem.mp.parameters.pl=[];
problem.mp.parameters.pu=[];
guess.mp.parameters=[];

% Bounds for linkage boundary constraints bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
problem.mp.constraints.bll.linear=[0 0];
problem.mp.constraints.blu.linear=[0 0];
problem.mp.constraints.blTol.linear=[0.01 0.01];

problem.mp.constraints.bll.nonlinear=[];
problem.mp.constraints.blu.nonlinear=[];
problem.mp.constraints.blTol.nonlinear=[];

% Get function handles
problem.mp.linkfunctions=@bclink;

% Store the necessary problem parameters used in the functions
problem.mp.data = [];

% Define different phases of OCP
[problem.phases{1},guess.phases{1}] = BangBangTwoPhase_Phase1(problem.mp, guess.mp);
[problem.phases{2},guess.phases{2}] = BangBangTwoPhase_Phase2(problem.mp, guess.mp);
phaseoptions{1}=problem.phases{1}.settings(20);
phaseoptions{2}=problem.phases{2}.settings(20);
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

blc_linear=[xf{1}(1)-x0{2}(1);xf{1}(2)-x0{2}(2)];
blc_nonlinear=[];
%------------- END OF CODE --------------
