function [ feasiSolve ] = checkConstraintReactivation( problem, solution, solution_old )
%checkConstraintReactivation - check if certain path constraints have become
%active again
%
% Syntax:  [ feasiSolve ] = checkConstraintReactivation( problem, solution, solution_old )
%
% Inputs:
%    problem - Defined in transcribeOCP
%    solution - the current solution
%    solution_old - the previous solution
%
% Outputs:
%    feasiSolve - boolean indicating whether the solve of feasiblity
%    problem is required
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

GlCheck=any(solution.ActiveConstraint.Glactive>solution_old.ActiveConstraint.Glactive);
GuCheck=any(solution.ActiveConstraint.Guactive>solution_old.ActiveConstraint.Guactive);

feasiSolve=any([GlCheck,GuCheck]);

end

