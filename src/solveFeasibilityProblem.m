function [ guess ] = solveFeasibilityProblem( problem_old, problem, options, data, guess)
%solveFeasibilityProblem - solve a feasiblity problem based on constraint violation, for use in external constraint handling 
%
% Syntax:  [ guess ] = solveFeasibilityProblem( problem_old, problem, options, data, guess)
%
% Inputs:
%    guess, options, data - Defined in transcribeOCP
%    problem - the modified problem 
%    problem_old - the original problem 
%
% Outputs:
%    guess - Modified initial guess, to ensure feasibility
%
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


guess_feas=guess;
n=length(problem.states.x0l);              % Number of states
m=length(problem.inputs.ul);               % Number of inputs

% Calculate the constraint violation and generate intial guess for the
% feasibility problem
    if strcmp(options.transcription,'globalLGR') || strcmp(options.transcription,'hpLGR')
        if isfield(options,'pdegree') 
            npd=options.pdegree*ones(1,options.nsegment);
            [npdu,~,npduidx]=unique(npd);
            nps=options.nsegment;
            M=sum(npd); 
        else
            npd=options.npsegment;
            [npdu,~,npduidx]=unique(npd);                  
            nps=length(options.npsegment);                
            M=sum(npd);                            
        end
        [ ~, tau_inc, ~, ~ ] = genTimeMeshLGR( options, data, nps, npdu, npd, npduidx, M );
        T=(guess.tf-data.t0)/2*tau_inc+(guess.tf+data.t0)/2;
        x_guess=zeros(M+1,n);u_guess=zeros(M,m);
        if isfield(guess,'TSeg_Bar')
            for i=1:n
                x_guess(:,i)=speval(guess.Xp,i,guess.TSeg_Bar,[T;guess.tf]);
            end
            for i=1:m
                u_guess(:,i)=speval(guess.Up,i,guess.TSeg_Bar,T);
            end
        else
            for i=1:n
                x_guess(:,i)=speval(guess.Xp,i,[T;guess.tf]);
            end
            for i=1:m
                u_guess(:,i)=speval(guess.Up,i,T);
            end
        end
        [ ConstraintError ] = max(calcConstraintViolation( data, guess, x_guess, [u_guess;u_guess(end,:)], problem, [T;guess.tf]));
        guess_feas.states=x_guess;
        guess_feas.inputs=[u_guess;u_guess(end,:)];
        guess_feas.time=[T;guess.tf];
    else
        if length(options.tau)~=1
            M=length(options.tau)+1;
        else    
            M=options.nodes;  
        end
        N=problem.inputs.N; 
        if N==0;N=M;end 
        if (strcmp(options.transcription,'hermite'))
            M=2*M-1;N=N*2-1;ns=2;
        else
            ns=1;
        end
        [ data, tau ] = genTimeMesh( options, data, ns, M );
        T=(guess.tf-data.t0)*[0;cumsum(tau)]*data.Nm/ns+data.k0;
        x_guess=zeros(M,n);u_guess=zeros(N,m);
        if isfield(guess,'TSeg_Bar')
            for i=1:n
                x_guess(:,i)=speval(guess.Xp,i,guess.TSeg_Bar,T);
            end
            for i=1:m
                u_guess(:,i)=speval(guess.Up,i,guess.TSeg_Bar,T);
            end
        else
            for i=1:n
                x_guess(:,i)=speval(guess.Xp,i,T);
            end
            for i=1:m
                u_guess(:,i)=speval(guess.Up,i,T);
            end
        end
        [ ConstraintError ] = max(calcConstraintViolation( data, guess, x_guess, u_guess, problem, T));
        guess_feas.states=x_guess;
        guess_feas.inputs=u_guess;
        guess_feas.time=T;
    end
    
    % Formulate the feasibility problem
    problem_feas=problem;
    problem_feas.data.mode.currentMode='Feasibility';
    ng=length(problem_feas.constraints.gl);
    problem_feas.data.mode.ng=ng;

    ConstraintError=ConstraintError(:);
    problem_feas.data.mode.np=length(ConstraintError)/2;
    problem_feas.time.tf_max=guess_feas.tf*1.05; %guess_feas.tf*1.05; 
    problem_feas.parameters.pl=zeros(1,length(ConstraintError));
    problem_feas.parameters.pu=inf*ones(1,length(ConstraintError));
    guess_feas.parameters=ConstraintError';

    options.start='Warm';
    problem_feas.constraints.gl=[-inf*ones(1,ng) problem_feas.constraints.gl];
    problem_feas.constraints.gu=[problem_feas.constraints.gu inf*ones(1,ng)];
    problem_feas.constraints.gTol=[problem_feas.constraints.gTol problem_feas.constraints.gTol];
    problem_feas.constraints.gActiveTime=[problem_feas.constraints.gActiveTime problem_feas.constraints.gActiveTime];
    
    
    options.ipopt.hessian_approximation='limited-memory';
    [infoNLP,data,options]=transcribeOCP(problem_feas,guess_feas,options); % Format for NLP solver
    [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
    [solution]=output(problem_old,solution,options,data,4);         % Output solutions
    
    %Use solution of feasibility as initial guess for next solve
    if status.status==0
        guess.states=solution.X;
        guess.Xp=solution.Xp;
        guess.inputs=solution.U;
        guess.Up=solution.Up;
        if isfield(solution,'TSeg_Bar')
            guess.TSeg_Bar=solution.TSeg_Bar;
        end
        guess.time=solution.T;
        guess.tf=solution.tf;
        guess.multipliers.lambda=solution.multipliers.lambda_1toN;
        nb=deal(data.sizes{6});
        if ng
            if strcmp(data.options.transcription,'hermite') || strcmp(data.options.transcription,'AutoDirect')
                guess.timeFull=solution.org.T;
                guess.multipliers.lambda_g=solution.multipliers.lambda_g(:,1:ng)+solution.multipliers.lambda_g(:,ng+1:end);
            else
                guess.timeFull=solution.T;
                guess.multipliers.lambda_g=solution.multipliers.lambda_g(:,1:ng)+solution.multipliers.lambda_g(:,ng+1:end);
            end
        end
        if nb
            guess.multipliers.lambda_b=solution.multipliers.lambda_b;
        end
    end
    
end

