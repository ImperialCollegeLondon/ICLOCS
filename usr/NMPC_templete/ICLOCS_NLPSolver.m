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
[problem,guess]=...;
    
% Updeate solution when time stamp reached
if simtime~=0 && simtime>solution_new.t_ref
    solution=solution_new;
end

re_opt=0;

% Conditions for forced re-optimization
if any(any(x0'>(problem.states.xu-problem.states.xConstraintTol))) %state close to upper bounds
    constVio=x0'-problem.states.xu;
    constVio(constVio<(-problem.states.xConstraintTol))=0;
    if any(constVio)
        re_opt=1;
    end
end
if any(any(x0'<(problem.states.xl+problem.states.xConstraintTol))) %state close to lower bounds
    constVio=problem.states.xl-x0';
    constVio(constVio<(-problem.states.xConstraintTol))=0;
    if any(constVio)
        re_opt=1;
    end
end
%add more conditions as needed (for event triggered reoptimization)

% Conditions for forced non-re-optimization
if simtime>solution.t_ref && simtime<(solution.t_ref+solution.computation_time)
    re_opt=0;
end
if solution.tf<opt_t_min
    re_opt=0;
end

if (simtime~=0 && ~re_opt) || mod(simtime,tstep)~=0 % For time instances without re-optimization
        crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
        crt_input=[speval(solution.Up,1,crt_time);speval(solution.Up,2,crt_time);speval(solution.Up,3,crt_time) ..] % for mutiple inputs, change as needed
        tfpu=[solution.tf-crt_time;0;crt_input] ; % Return the interpolated solution
else
    options= settings_h(N);                  % Get options and solver settings 
    if simtime~=0 % For later instances
        % If appromated value exceed the state bounds, set the value to
        % state bound limits if within user-defined tolerance, otherwise
        % terminate with an error
        opt_step=solution.computation_time;
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
        idx_new=solution.T>(opt_step);
        T_new=[0;solution.T(idx_new)-(opt_step)];
        guess.tf=solution.tf-(opt_step);
        guess.time=T_new;
        guess.states=[x0';solution.X(idx_new,:)];
        crt_input=[speval(solution.Up,1,opt_step) speval(solution.Up,2,opt_step) speval(solution.Up,3,opt_step) .. ] % for mutiple inputs, change as needed
        guess.inputs=[crt_input;solution.U(idx_new,:)];
    
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
        solution_new.computation_time=computation_time; %set (next) update interval to the time needed for this compuation
        solution_new.t_ref=simtime; % Give the solution a time stamp
        % still return interpolated results from previous solution (as the new one not yet available due to comp. delay)
        crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
        crt_input=[speval(solution.Up,1,crt_time);speval(solution.Up,2,crt_time);speval(solution.Up,3,crt_time) ..] % for mutiple inputs, change as needed
        tfpu=[solution.tf-opt_step;computation_time;crt_input] ; % Return the interpolated solution
    else
        % Assume solution of the first solve already available at t=0
        solution_new.computation_time=0;
        solution_new.t_ref=0;
        tfpu=[solution_new.tf;solution_new.computation_time;solution_new.U(1,:)'];
        solution=solution_new;
    end

end

end



