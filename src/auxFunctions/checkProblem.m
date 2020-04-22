function checkProblem(problem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
