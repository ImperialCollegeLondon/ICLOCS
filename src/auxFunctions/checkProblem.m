function checkProblem(problem)
%checkProblem Validate key dimensions in an ICLOCS problem structure.
%
% Syntax:
%   checkProblem(problem)
%
% Input:
%   problem - ICLOCS problem structure returned by a problem-definition file.
%
% Output:
%   None.
%
% Notes:
%   Throws an error when state, input, path-constraint, or boundary-constraint
%   bounds and tolerances have inconsistent dimensions. Also checks for the
%   ICLOCS 2.5 problem.data.ng_eq convention.

if length(problem.states.x0l)==length(problem.states.x0u) && length(problem.states.x0l)==length(problem.states.xl) && length(problem.states.x0l)==length(problem.states.xu) && length(problem.states.x0l)==length(problem.states.xErrorTol_local) && length(problem.states.x0l)==length(problem.states.xErrorTol_integral) && length(problem.states.x0l)==length(problem.states.xConstraintTol) && length(problem.states.x0l)==length(problem.states.xfl) && length(problem.states.x0l)==length(problem.states.xfu)
else
    error('Please ensure all information in problem formulation regarding the state variables and their error tolerances have been correctly configured')
end

if length(problem.inputs.ul)==length(problem.inputs.uu) && length(problem.inputs.ul)==length(problem.inputs.uConstraintTol)
else
    error('Please ensure all information in problem formulation regarding the input variables and their error tolerances have been correctly configured')
end

if problem.constraints.ng_eq==length(problem.constraints.gTol_eq)
else
    error('Please ensure all information in problem formulation regarding the equality path constraints and their error tolerances have been correctly configured')
end

if length(problem.constraints.gl)==length(problem.constraints.gu) && length(problem.constraints.gl)==length(problem.constraints.gTol_neq)
else
    error('Please ensure all information in problem formulation regarding the inequality path constraints and their error tolerances have been correctly configured')
end

if length(problem.constraints.bl)==length(problem.constraints.bu) && length(problem.constraints.bl)==length(problem.constraints.bTol)
else
    error('Please ensure all information in problem formulation regarding the boundary constraints and their error tolerances have been correctly configured')
end

if ~isfield(problem.data,'ng_eq')
    error('Please follow the instructions on https://github.com/ImperialCollegeLondon/ICLOCS to update your problem to run with ICLOCS version 2.5')
end
