function [ problem,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, iter )
%selectAppliedConstraint - Select the constraints to be applied, for external constraint handling 
%
% Syntax:  [ problem,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, iter )
%
% Inputs:
%    problem, guess, options, data - Defined in transcribeOCP
%    solutionHistory - Structure containing all past solutions
%    iter - current re-solve iteration
%
% Outputs:
%    problem,guess,options - Modified problem to be solved
%
% Other m-files required: identifyConstActive, checkConstraintReactivation, solveFeasibilityProblem
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

%------------- BEGIN CODE --------------


% pull the previous solution
solution=solutionHistory{iter-1};

if (strcmp(options.discretization,'globalLGR')) || (strcmp(options.discretization,'hpLGR'))
    ng_eq=deal(data.sizes{19});
else
    ng_eq=deal(data.sizes{15});
end

ng=deal(data.sizes{5});

%% Remove variable simple bounds if never become active
ActiveConstraint=solution.ActiveConstraint;
problem_old = problem;

%State bounds. xl=< x <=xu
problem.states.xl(ActiveConstraint.Xlactive)=problem_old.states.xl(ActiveConstraint.Xlactive);
problem.states.xl(~ActiveConstraint.Xlactive)=-inf;
problem.states.xu(ActiveConstraint.Xuactive)=problem_old.states.xu(ActiveConstraint.Xuactive);
problem.states.xu(~ActiveConstraint.Xuactive)=inf;

% State rate bounds. xrl=< x <=xru
if isfield(solution,'dX') && isfield(problem.states,'xrl')
    problem.states.xrl(ActiveConstraint.Xrlactive)=problem_old.states.xrl(ActiveConstraint.Xrlactive);
    problem.states.xrl(~ActiveConstraint.Xrlactive)=-inf;
    problem.states.xru(ActiveConstraint.Xruactive)=problem_old.states.xru(ActiveConstraint.Xruactive);
    problem.states.xru(~ActiveConstraint.Xruactive)=inf;
end


%% Remove potentially inactive path constraints
if ng
    % Construct multiplier infomration for the original problem
    guess_lambda_g=zeros(size(solution.multipliers.lambda_g,1),length(ActiveConstraint.Glactive));
    if isfield(guess.multipliers,'lambda_g') && isfield(data.data,'gFilter')
        guess_lambda_g(:,~data.data.gFilter)=guess.multipliers.lambda_g;
    else
        guess_lambda_g=guess.multipliers.lambda_g;
    end

    % For all path constraints, identify at which time intervals this
    % will happen
    idxGActive0=true(1,size(guess_lambda_g,2));
    
    for i=ng_eq+1:size(guess_lambda_g,2)
        gActiveTime{i} = identifyConstActive( problem_old, options, solution, guess_lambda_g(:,i), ng, i );
        if isempty(gActiveTime{i}) % No time intervals can be indentified for constraint being active
            idxGActive0(i)=false; % consider that constraint potentially inactive (criterion 1)
        end
    end
    
    % Path constraint considered potentially active if either the upper or
    % the lower bounds has been reached (criterion 2)
    idxGActive1=ActiveConstraint.Glactive | ActiveConstraint.Guactive;

    % Path constraint considered potentially active if singificant
    % variations in the costate values are observed (criterion 3)
    multiG_std= abs(std(guess_lambda_g(2:end-1,ng_eq+1:end)));
    idxGActive2=[true(1,ng_eq),(multiG_std)> mean(multiG_std)];

    % Path constraint considered active if criterion 2 is fulfilled, or
    % criterion 1 and 3 are both fulfilled. idxjoint containts a list of
    % potentially inactive constraints.
    idxjoint = idxGActive2 & idxGActive0;
    idxjoint=~(idxGActive1 | idxjoint);

    % Remove potentially inactive constraints from the OCP
    problem.constraints.gl(idxjoint(ng_eq+1:end))=[];
    problem.constraints.gu(idxjoint(ng_eq+1:end))=[];
    problem.constraints.gTol_neq(idxjoint(ng_eq+1:end))=[];

    guess.multipliers.lambda_g=guess_lambda_g(:,~idxjoint);

    gActiveTime(:,idxjoint)=[];
    problem.constraints.g_neq_ActiveTime=gActiveTime(:,ng_eq+1:end);
    
    if isfield(options.ECH,'noConstIntv') && options.ECH.noConstIntv
        problem.constraints=rmfield(problem.constraints,'g_neq_ActiveTime');
    end

    problem.data.gFilter=idxjoint;
    
end

problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

if isfield(data.data,'Xscale')
    problem.states.scales=data.data.Xscale;
    problem.states.scales_back=data.data.Xscale_back;
    problem.states.shifts=data.data.Xshift;
    problem.inputs.scales=data.data.Uscale;
    problem.inputs.scales_back=data.data.Uscale_back;
    problem.inputs.shifts=data.data.Ushift;
end

% For later iterations, check if constraints removed earlier has become
% potentially active again. If ture, solve the feasibility problem first
if iter>2
    solution_old=solutionHistory{iter-2};
    [ feasiSolve ] = checkConstraintReactivation( problem, solution, solution_old );
    if feasiSolve
        [ guess ] = solveFeasibilityProblem( problem_old, problem, options, data, guess);
    end
end

end

