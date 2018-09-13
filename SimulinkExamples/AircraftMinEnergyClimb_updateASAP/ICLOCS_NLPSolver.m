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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

global infoNLP data solution N opt_t_min tstep
persistent solution_new

simtime=tpx0(1); %Current simulation time (global)
x0=tpx0(2:end); %Current state

% Fetch the problem definition
[problem,guess]=MinEnergyClimbBryson;          

% Updeate solution when time stamp reached
if simtime~=0 && simtime>solution_new.t_ref
    solution=solution_new;
end

re_opt=0;

% Conditions for forced re-optimization
if any(any(x0'>(problem.states.xu-problem.states.xConstraintTol)))
    constVio=x0'-problem.states.xu;
    constVio(constVio<(-problem.states.xConstraintTol))=0;
    if any(constVio)
        re_opt=1;
    end
end
if any(any(x0'<(problem.states.xl+problem.states.xConstraintTol)))
    constVio=problem.states.xl-x0';
    constVio(constVio<(-problem.states.xConstraintTol))=0;
    if any(constVio)
        re_opt=1;
    end
end

% Conditions for forced non-re-optimization
if simtime>solution.t_ref && simtime<(solution.t_ref+solution.computation_time)
    re_opt=0;
end
if solution.tf<opt_t_min
    re_opt=0;
end

if (simtime~=0 && ~re_opt) || mod(simtime,tstep)~=0 % For time instances without re-optimization
    crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
    tfpu=[solution.tf-crt_time;0;speval(solution.Up,1,crt_time)] ; % Return the interpolated solution
else
    options= settings_h(N);                  % Get options and solver settings 
    if simtime~=0 % For later instances
        
        opt_step=solution.computation_time;
        
        % If appromated value exceed the state bounds, set the value to
        % state bound limits if within user-defined tolerance, otherwise
        % terminate with an error
        if any(any(x0'>problem.states.xu))
            constVio=x0'-problem.states.xu;
            constVio(constVio<0)=0;
            if all(constVio<problem.states.xConstraintTol)
                x0(constVio>0)=problem.states.xu(constVio>0);
            else
                error('Simulation going out of bound');
            end
        end
        if any(any(x0'<problem.states.xl))
            constVio=problem.states.xl-x0';
            constVio(constVio<0)=0;
            if all(constVio<problem.states.xConstraintTol)
                x0(constVio>0)=problem.states.xl(constVio>0);
            else
                error('Simulation going out of bound');
            end
        end
        
        % Update the initial condition
        problem.states.x0=x0';
        problem.states.x0l=x0';
        problem.states.x0u=x0';
        
        % Update the initial guess
        if solution.tf > opt_step
            idx_new=solution.T>(opt_step);
            T_new=[0;solution.T(idx_new)-(opt_step)];
            guess.tf=solution.tf-(opt_step);
            guess.time=T_new;
            guess.states=[x0';solution.X(idx_new,:)];
            guess.inputs=[speval(solution.Up,1,opt_step);solution.U(idx_new,:)];
        else
            guess.states(1,:)=x0';
        end
        
        % Update the new mesh
        idx_new_2=solution.T>opt_step;
        T_new_2=[0;solution.T(idx_new_2)-opt_step];
        if solution.tf>opt_step && length(T_new_2)>=3
            options.tau=diff(T_new_2./guess.tf);
            options.nodes=length(options.tau)+1;
        else
            T_new_2=linspace(0,guess.tf,3);
            T_new_2(end)=guess.tf;
            options.tau=diff(T_new_2./guess.tf)';
            options.nodes=length(options.tau)+1;  
        end
    end

    % Solve the optimization problem
    maxAbsError=1e9;maxAbsConstraintError=1e9;
    i=1; imax=50;
    computation_time=0;
    while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax && solution.status.status~=2
        [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
        [solution_new,solution_new.status,data] = solveNLP(infoNLP,data);      % Solve the NLP
        [solution_new]=output(problem,solution_new,options,data,0);          % Output solutions
        computation_time=computation_time+solution_new.computation_time;
        maxAbsError=max(abs(solution_new.Error));
        maxAbsConstraintError=max(solution_new.ConstraintError);

        [ options, guess ] = doMeshRefinement( options, problem, guess, data, solution_new, i );
        i=i+1;
    end

    % Return the solution
    if simtime~=0
        solution_new.computation_time=computation_time;
        solution_new.t_ref=simtime;
        crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
        tfpu=[solution.tf-crt_time;computation_time;speval(solution.Up,1,crt_time)] ; % Return the interpolated solution
    else
        solution_new.computation_time=0;
        solution_new.t_ref=0;
        tfpu=[solution_new.tf;solution_new.computation_time;solution_new.U(1,:)'];
        solution=solution_new;
    end

end

end



