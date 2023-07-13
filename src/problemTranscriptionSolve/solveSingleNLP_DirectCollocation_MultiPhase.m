function [solution,status,data]=solveSingleNLP_DirectCollocation_MultiPhase(NLP,data)

%solveSingleNLP_DirectCollocation - Solve a single multi-phase nonlinear program and
%return the solution (direct collocation method)
%
% Syntax:  [solution,status,data]=solveSingleNLP_DirectCollocation(NLP,data)
%
% Inputs:
%    NLP     - Information required by the NLP solver
%    data    - Data passed to the functions evaluated during optimization
%
% Outputs:
%    solution - Structure containing the optimal final time, parameters,
%               states, controls and Lagrange multipliers.
%    status - variable containing the exit condition of the solver. See the
%             fmincon and ipopt description 
%
% Other m-files required: ipopt.mex and or fmincon.m
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

% Solve the nonlinear program:




  

switch(data.mpdata.options.mp.NLPsolver)

    case{'ipopt'} 
       if isfield(data.mpdata.multipliers,'lambda') && ~isempty(data.mpdata.multipliers.lambda)
            opt.lambda=data.mpdata.multipliers.lambda;              % Solve the NLP using IPOPT
            data.options.ipopt.warm_start_init_point ='yes'; 
       else
            data.mpdata.options.ipopt.warm_start_init_point ='no'; 
       end
       opt.ipopt=data.mpdata.options.mp.ipopt;
       opt.zl=NLP.mpinfoNLP.zl;
       opt.zu=NLP.mpinfoNLP.zu;
       opt.cl=NLP.mpinfoNLP.cl;                    % Lower bounds on constraints.
       opt.cu=NLP.mpinfoNLP.cu;                      % Upper bounds on constraints
       opt.lb=NLP.mpinfoNLP.zl;                       % Lower bound on the variables.
       opt.ub=NLP.mpinfoNLP.zu;                      % Upper bound on the variables.
       opt.auxdata=data;
       
       tA=tic;
       [z, info] = ipopt_auxdata(NLP.mpinfoNLP.z0,data.mpdata.funcs,opt);
       tB=toc(tA);

       solution.mp.z_org=z;
       solution.mp.multipliers.lambda=info.lambda;
       solution.mp.multipliers.zl=info.zl;
       solution.mp.multipliers.zu=info.zu;
       status=info;
       
       
%        if  ((strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))) && data.options.reorderLGR
%            solution.z_org=z;
%            z=z(data.reorder.z_idx_back);
%            info.zl=info.zl(data.reorder.z_idx_back);
%            info.zu=info.zu(data.reorder.z_idx_back);
%            info.lambda=info.lambda(data.reorder.vert_idx_back);
%        end
%        
%        solution.multipliers.zl=info.zl;
%        solution.multipliers.zu=info.zu;
%        solution.multipliers.lambda=info.lambda;
%        status=info;
       
    
    case{'fmincon'}                               % Solve the NLP using FMINCON
        error('use of fmincon not yet supported for multi-phase problem, please use IPOPT')
        prompt = 'WARNING: Use of fmincon is NOT ADVISABLE because the solution speed can be significantly slower with a higher chance of failure. \n We would recommend to use fmincon only for sanity checks. \n For solving, please consider use IPOPT (see http://www.ee.ic.ac.uk/ICLOCS/Downloads.html). \n Press any key to continue...\n ';
        input(prompt)
        tA=tic;
        [z,cost,status,output,multipliers]=fmincon(@(z)fminCost(z,data),NLP.z0,...
            [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.fmincon);
        tB=toc(tA);
        if  ((strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))) && data.options.reorderLGR
           solution.z_org=z;
           z=z(data.reorder.z_idx_back);
           multipliers.eqnonlin=multipliers.eqnonlin(data.reorder.vert_idx_back);
        end
        ni=output.iterations;
        solution.multipliers.lambda=multipliers.eqnonlin;
        solution.iterates=ni;
        solution.cost=cost;

    case{'NOMAD'}                               % Solve the NLP using FMINCON
        tA=tic;
        
        obj_type = cell(length(NLP.cu)*2,1);
        obj_type(:) = {'PB'};
        if isfield(data.options,'NOMAD')
           opts = nomadset('display_degree',2,'bb_output_type',obj_type,'initial_mesh_size','50000','h_min',1e-03,'f_target',data.options.NOMAD.f_target);
        else
           opts = nomadset('display_degree',2,'bb_output_type',obj_type,'initial_mesh_size','50000','h_min',1e-03,'h_norm','L2');
        end
        xtype    = cell(length(NLP.z0),1);
        xtype(:) = {'C'};
        xtype=char(xtype);
        [z,cost,status,output.iter,output.nfval] = nomad(@(z)MADSCost(z,data),NLP.z0,NLP.zl,NLP.zu,@(z)MADSConst(z,data,NLP),[NLP.cu;-NLP.cl],xtype,opts);
%         [z,cost,status,output.iter,output.nfval] = nomad(@(z)MADSCost(z,data),NLP.z0,NLP.zl,NLP.zu,[],[],xtype,opts);

%         [z,cost,status,output]=patternsearch(@(z)fminCost(z,data),NLP.z0,...
%             [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.MADS);
        tB=toc(tA);
        if  ((strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))) && data.options.reorderLGR
           solution.z_org=z;
           z=z(data.reorder.z_idx_back);
        end
        ni=output.iter;
        solution.iterates=ni;
        solution.cost=cost;
        
    case{'builtinSQP'}
        SQPoptions=data.options.SQP;
        if strcmp(data.options.ipopt.hessian_approximation,'exact') 
            SQPoptions.realHess = 1;
        else
            SQPoptions.realHess = 0;
        end
        
        if strcmp(SQPoptions.QPSolver,'quadprog') 
            qp_options = optimoptions(@quadprog,'disp','off');
            qpSolver = @(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)qps_ip(H,f,A,b,Aeq,beq,lb,ub,x0,qp_options,Ac);
        elseif strcmp(SQPoptions.QPSolver,'qpOASES') 
            qpSolver = @(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)qps_qpOASES(H,f,A,b,Aeq,beq,lb,ub,x0,Ac);
        else
            error('invalid choice of QP solver for the built-in SQP');
        end
        
        bounds.xl=NLP.zl;
        bounds.xu=NLP.zu;
        bounds.cl=NLP.cl;                    % Lower bounds on constraints.
        bounds.cu=NLP.cu;                      % Upper bounds on constraints
        bounds.bl=[];                       
        bounds.bu=[];                      
        
        tA=tic;
        [z, info] = LISSQP(data.funcs,[],bounds, qpSolver, SQPoptions, data, NLP.z0);
        tB=toc(tA);
        
        if  ((strcmp(data.options.discretization,'globalLGR')) || (strcmp(data.options.discretization,'hpLGR'))) && data.options.reorderLGR
           solution.z_org=z;
           z=z(data.reorder.z_idx_back);
           info.lambda=info.lambda(data.reorder.vert_idx_back);
        end
        
        solution.multipliers.lambda=info.lambda;
        status=[];
    case{'worhp'}
        
        tA=tic;
        
        dataworhp = worhp('worhp.xml');
        z=NLP.z0;

        La=ones(size(z));
        if ~isempty(data.multipliers)
            Mu=-data.multipliers.lambda';
        else
            Mu=ones(size(NLP.cu));
        end
        [actions, scaleF, status] = worhp(dataworhp, NLP.zl, z, NLP.zu, NLP.cl, zeros(size(NLP.cu)), NLP.cu, La, Mu, 1:length(z), data.jacStruct_GRow, data.jacStruct_GCol, data.jacStruct_HMRow, data.jacStruct_HMCol);
        
        if strcmp(data.options.discretization,'globalLGR') || strcmp(data.options.discretization,'hpLGR')
            while (actions.continue)

               if (actions.worhp)
                 [actions, scaleF, status, zNew, LaNew, MuNew] = worhp(dataworhp, z, La, Mu);
                 % DO NOT REMOVE THESE LINES
                 % Use updated variables from WORHP from now on
                 z  = zNew;
                 La = LaNew;
                 Mu = MuNew;
                 % DO NOT REMOVE THESE LINES
                 data.lambda=[-Mu'];data.sigma=1;
               end

                if (actions.output)
                  [actions, scaleF, status] = worhp(dataworhp);
                end

                if (actions.f)
                    F = scaleF*directCollocationLGR('cost',z',data);
                    [actions, scaleF, status] = worhp(dataworhp, F, 'evalF');
                end

                if (actions.g)
                    G = directCollocationLGR('const',z',data)';
                    [actions, scaleF, status] = worhp(dataworhp, G, 'evalG');
                end

                if (actions.df)
                    DFVAL = scaleF * full(directCollocationLGR('gradCost',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, DFVAL, 'evalDF');
                end

                if (actions.dg)
                    DGVAL = full(directCollocationLGR('jacConst',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, DGVAL, 'evalDG');
                end

                if (actions.hm)
                    HMVAL = full(directCollocationLGR('hessian',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, HMVAL, 'evalHM');
                end

                if (actions.fidif)
                  [actions, scaleF, status] = worhp(dataworhp, z, La, F, G, Mu);
                end
            end
        else
            while (actions.continue)

               if (actions.worhp)
                 [actions, scaleF, status, zNew, LaNew, MuNew] = worhp(dataworhp, z, La, Mu);
                 % DO NOT REMOVE THESE LINES
                 % Use updated variables from WORHP from now on
                 z  = zNew;
                 La = LaNew;
                 Mu = MuNew;
                 % DO NOT REMOVE THESE LINES
                 data.lambda=[-Mu'];data.sigma=1;
               end

                if (actions.output)
                  [actions, scaleF, status] = worhp(dataworhp);
                end

                if (actions.f)
                    F = scaleF*directCollocation('cost',z',data);
                    [actions, scaleF, status] = worhp(dataworhp, F, 'evalF');
                end

                if (actions.g)
                    G = directCollocation('const',z',data)';
                    [actions, scaleF, status] = worhp(dataworhp, G, 'evalG');
                end

                if (actions.df)
                    DFVAL = scaleF * full(directCollocation('gradCost',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, DFVAL, 'evalDF');
                end

                if (actions.dg)
                    DGVAL = full(directCollocation('jacConst',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, DGVAL, 'evalDG');
                end

                if (actions.hm)
                    HMVAL = full(directCollocation('hessian',z',data));
                    [actions, scaleF, status] = worhp(dataworhp, HMVAL, 'evalHM');
                end

                if (actions.fidif)
                  [actions, scaleF, status] = worhp(dataworhp, z, La, F, G, Mu);
                end
            end
        end
        
        % Call Worhp again for final output
        worhp(dataworhp);
        z=z';
        tB=toc(tA);
        solution.z_org=z;
        solution.multipliers.lambda=-Mu';

    otherwise
        disp('Unknown NLP solver. Check spelling.');

end


% Store the results in solution structure:
%solution.status=status;
solution.mp.computation_time=tB;

try
    for i=1:length(data.phasedata)
        z_phase=zeros(data.phasedata{i}.nz,1);
        sol.multipliers.lambda=zeros(data.phasedata{i}.nConst,1);
        z_phase(data.phasedata{i}.zidx.org.z,1)=z(data.phasedata{i}.zidx.mp.z,1);
        sol.multipliers.lambda(data.phasedata{i}.zidx.org.const,1)=solution.mp.multipliers.lambda(data.phasedata{i}.zidx.mp.const,1);
        [solution.phaseSol{i}] = postProcessSolution(z_phase, data.phasedata{i}, sol);
    end
catch
    error('Error encountered when post-processing the solution. Please ensure the NLP solve has been terminated successfully, and the error tolerances have been correctly configured');
end



%------------- END OF CODE --------------


end

%% Functions needed for the built-in SQP
function [x, lambdaQP, flag,Ac] = qps_ip(H,f,A,b,Aeq,beq,lb,ub,x0,opt,Ac)
[x,t1,flag,t2,l_QuadProg] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opt);
lambdaQP = [ l_QuadProg.eqlin; l_QuadProg.ineqlin; l_QuadProg.lower; l_QuadProg.upper];
Ac = [];
end

function [x, lambdaQP, flag,Ac] = qps_qpOASES(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)

At  = [Aeq; A];
ubA = [beq ; b];
lbA = [beq ; -1*inf*ones(size(b))];

m1 = size(Aeq,1);
m2 = size(A,1);
m3 = length(lb);

WorkingSetGuess = qpOASES_auxInput('guessedWorkingSetB',Ac(1:length(x0)),'guessedWorkingSetC',Ac(length(x0)+1:end));
[x,fval,flag,iter,lambdaQP,auxOut] = qpOASES(H,f,At,lb,ub,lbA,ubA,qpOASES_options,WorkingSetGuess);

% -1 indicates lower bound is active
% 0 indicates that no bound is active
% 1 indicates that the upper bound is active 
W1 = auxOut.workingSetB;
W2 = auxOut.workingSetC;

% Generate masks for lagrange multipliers from working set
W1_l = -1*W1;
W1_l(W1_l<0) = 0;

W1_u= W1;
W1_u(W1_u<0) = 0;

% Get multipliers and apply active set mask for inequalities
lambdaQP = [-1*lambdaQP(m3+1:m3+m1); ...
    -1*lambdaQP(m3+m1+1 : end) .* (W2(m1+1:m1+m2)); ...
    [lambdaQP(1:m3) ; lambdaQP(1:m3)].*[W1_l ; -1*W1_u]];
Ac = [W1;W2];

end