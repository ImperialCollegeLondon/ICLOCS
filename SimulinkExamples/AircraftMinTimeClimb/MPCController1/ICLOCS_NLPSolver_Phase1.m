function [ tfpu ] = ICLOCS_NLPSolver_Phase1(tpx0)

%ICLOCS_NLPSolver - Implement optimal control problem with ICLOCS in
%closed-loop with Simulink
%
% Syntax:  [ tfpu ] = ICLOCS_NLPSolver_Phase1(tpx0)
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

global solution N_node opt_t_min tstep settings_func opt_step
persistent solution_new

simtime=tpx0(1); %Current simulation time (global)
x0=tpx0(2:end); %Current state

% Fetch the problem definition
[problem,guess]=MinTimeClimbBryson_Phase1;          

% Update solution when time stamp reached
if simtime~=0 && simtime>=(solution_new.t_ref+solution_new.computation_time)
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
if simtime==0 || simtime>(solution.t_ref+solution.tf)
    re_opt=1;
end

% Standard Update
if simtime~=0 && simtime>=solution.t_ref+solution.computation_time
    re_opt=1;
end

% Event Based Update
% if simtime>(solution.t_ref+20)
%     re_opt=1;
% end


% Conditions for forced non-re-optimization
if simtime>solution.t_ref && simtime<(solution_new.t_ref+solution_new.computation_time)
    re_opt=0;
end
if solution.tf<opt_t_min && simtime<(solution.t_ref+solution.tf)
    re_opt=0;
end

if (simtime~=0 && ~re_opt) || mod(simtime,tstep)~=0 % For time instances without re-optimization
    crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
    tfpu=[solution.tf-crt_time;0;speval(solution,'U',1,crt_time)] ; % Return the interpolated solution
else
    options= settings_func(N_node);                  % Get options and solver settings 
    if simtime~=0 % For later instances
        
        opt_st=solution.computation_time;
        
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
        if solution.tf > opt_st
            idx_new=solution.T>(opt_st);
            T_new=[0;solution.T(idx_new)-(opt_st)];
            guess.tf=solution.tf-(opt_st);
            guess.time=T_new;
            guess.states=[x0';solution.X(idx_new,:)];
            guess.inputs=[speval(solution,'U',1,opt_st);solution.U(idx_new,:)];
        else
            guess.states(1,:)=x0';
        end
        
        if problem.time.tf_max-simtime>0.1
            problem.time.tf_max=problem.time.tf_max-simtime;
        else
            problem.time.tf_max=0.1;
        end
       
        guess.tf(guess.tf<problem.time.tf_min)=problem.time.tf_min;
        guess.tf(guess.tf>problem.time.tf_max)=problem.time.tf_max;
        
        % Update the new mesh
        idx_new_2=solution.T>opt_st;
        T_new_2=[0;solution.T(idx_new_2)-opt_st];
        if solution.tf>opt_st && length(T_new_2)>=10
            options.tau=diff(T_new_2./T_new_2(end));
            options.nodes=length(options.tau)+1;
        else
            T_new_2=linspace(0,guess.tf,10);
            T_new_2(end)=guess.tf;
            options.tau=diff(T_new_2./guess.tf)';
            options.nodes=length(options.tau)+1;  
        end
    end

    % Solve the optimization problem
    [solution_new,MRHistory]=solveMyProblem( problem,guess,options);
    computation_time=sum(MRHistory.timeHistory);
    
    % Return the solution
    if simtime~=0
        if opt_step==0
            solution_new.computation_time=computation_time; %Update at every time step
        else
            solution_new.computation_time=opt_step; %Update at every time step
        end
        solution_new.t_ref=simtime;
        crt_time=simtime-solution.t_ref; % Current time in the frame of previous solution
        tfpu=[solution.tf-crt_time;computation_time;speval(solution,'U',1,crt_time)] ; % Return the interpolated solution
    else
        solution_new.computation_time=opt_step;
        solution_new.t_ref=0;
        tfpu=[solution_new.tf;solution_new.computation_time;solution_new.U(1,:)'];
        solution=solution_new;
    end

end

end



