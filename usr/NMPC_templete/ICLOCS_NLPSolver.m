function [ tfpu ] = ICLOCS_NLPSolver(tpx0)

%ICLOCS_NLPSolver - Implement optimal control problem with ICLOCS in
%closed-loop with Simulink
%
% Syntax:  [ tfpu ] = ICLOCS_NLPSolver(tpx0)
%
% Input:
%    tpx0 - vector in format [current time; current state]
%
% Output:
%    tfpu - vector in format [terminal time; computation time; input to implement]
%
% Other m-files required: none
% Subfunctions: none
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

x0=tpx0(2:end); %Current state
global opt_step opt_t_min

if  ConditionForController2
    opt_step=...; % time interval for re-optimization for controller 2, set to zero for update ASAP
    opt_t_min=0; % stop re-optimizing if solution.tf<opt_t_min
    [ tfpu ] = ICLOCS_NLPSolver_Phase2(tpx0);
...
else % default controller
    opt_step=0.01; % time interval for re-optimization, set to zero for update ASAP
    opt_t_min=0; % stop re-optimizing if solution.tf<opt_t_min
    [ tfpu ] = ICLOCS_NLPSolver_Phase1(tpx0);
end

end



