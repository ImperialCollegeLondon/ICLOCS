function [solution,status,data]=solveNLP(NLP,data)

%SOLVENLP - Solve the nonlinear program and return the solution
%
% Syntax:  solution = solveNLP(NLP,options)
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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

% Solve the nonlinear program:



if strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce]=deal(data.sizes{:});
    if data.options.reorderLGR
        NLP.z0(1:n)=data.x0t;
    else
        NLP.z0(1:M+1:M*n)=data.x0t;
    end
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce]=deal(data.sizes{:});
    NLP.z0(nt+np+1:n+nt+np)=data.x0t;
end
nrc=nrcl+nrcu+nrce;

  

switch(data.options.NLPsolver)

    case{'ipopt'} 
       opt=data.multipliers;              % Solve the NLP using IPOPT
       opt.ipopt=data.options.ipopt;
       opt.zl=NLP.zl;
       opt.zu=NLP.zu;
       opt.cl=NLP.cl;                    % Lower bounds on constraints.
       opt.cu=NLP.cu;                      % Upper bounds on constraints
       opt.lb=NLP.zl;                       % Lower bound on the variables.
       opt.ub=NLP.zu;                      % Upper bound on the variables.
       opt.auxdata=data;
       
       tA=tic;
       [z, info] = ipopt_auxdata(NLP.z0,data.funcs,opt);
       tB=toc(tA);
       
       if  ((strcmp(data.options.transcription,'globalLGR')) || (strcmp(data.options.transcription,'hpLGR'))) && data.options.reorderLGR
           solution.z_org=z;
           z=z(data.reorder.z_idx_back);
           info.zl=info.zl(data.reorder.z_idx_back);
           info.zu=info.zu(data.reorder.z_idx_back);
           info.lambda=info.lambda(data.reorder.vert_idx_back);
       end
       
       solution.multipliers.zl=info.zl;
       solution.multipliers.zu=info.zu;
       solution.multipliers.lambda=info.lambda;
       status=info.status;
       
    
    case{'fmincon'}                               % Solve the NLP using FMINCON
        
        tA=tic;
        [z,cost,status,output,multipliers]=fmincon(@(z)fminCost(z,data),NLP.z0,...
            [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.fmincon);
        tB=toc(tA);
        if  ((strcmp(data.options.transcription,'globalLGR')) || (strcmp(data.options.transcription,'hpLGR'))) && data.options.reorderLGR
           solution.z_org=z;
           z=z(data.reorder.z_idx_back);
           multipliers.eqnonlin=multipliers.eqnonlin(data.reorder.vert_idx_back);
        end
        ni=output.iterations;
        solution.multipliers.lambda=multipliers.eqnonlin;
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
        
        if  ((strcmp(data.options.transcription,'globalLGR')) || (strcmp(data.options.transcription,'hpLGR'))) && data.options.reorderLGR
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
        
        if strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')
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
solution.computation_time=tB;

% solution.z=NLP.z0;
solution.z=z;


if strcmp(data.options.transcription,'globalLGR') || strcmp(data.options.transcription,'hpLGR')
    if nt;solution.tf=z(end); else solution.tf=data.tf; end
    if np;solution.p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np); else solution.p=[]; end
    t0=z(end-nt+1);
    solution.t0=t0;

    solution.X=reshape(data.map.Vx*z,M+1,n);
    solution.x0=solution.X(1,:)';
    solution.U=reshape(data.map.Vu*z,M,m);
    
    if data.options.adaptseg==1
    t_segment=z(end-nt+1:end)';
    
    tau_segment = ((1-(-1))*t_segment-(t_segment(1)+t_segment(end))) /(t_segment(end)-t_segment(1));
    
    tau_inc=[];

    tau_seg=cell(1,length(tau_segment)-1);
    for i=1:length(tau_segment)-1
        tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [data.map.LGR.points{npduidx(i)}; 1]) / (1-(-1));
        tau_inc=[tau_inc; tau_seg{i}(1:end-1)];
    end
        
    tau=diff(tau_inc);
    data.tau_inc=tau_inc;
    data.tau_segment=tau_segment;
    data.tau=tau;
    end
    
    solution.T=(solution.tf-t0)/2*data.tau_inc+(solution.tf+t0)/2; 
    solution.multipliers.lambdaNLP=solution.multipliers.lambda;
    solution.multipliers.lambda_1toN=reshape(solution.multipliers.lambda(1:M*n),M,n);
    D_N=diag(1./data.map.w);
    solution.multipliers.lambda=D_N*solution.multipliers.lambda_1toN;
    D_Np1=data.map.LGR.diff_matrix{end}(:,end);
    solution.multipliers.lambda(end+1,:)=D_Np1'*solution.multipliers.lambda_1toN(end-length(data.map.LGR.points{end})+1:end,:);
    if ng
        solution.multipliers.lambda_g=reshape(solution.multipliers.lambdaNLP((M*n+1):M*(n+ng)),M,ng);
    end
    if nrc
        solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP((M*(n+ng)+1):M*(n+ng)+nrc);
    end
    if nb
        solution.multipliers.lambda_b=solution.multipliers.lambdaNLP((M*(n+ng)+nrc+1):end);
    end

else
    
    if nt
        solution.tf=z(1);
    else
        solution.tf=data.tf; 
    end
    
    if np
        solution.p=z(nt+1:nt+np);
    else
        solution.p=[]; 
    end
    t0=data.t0;

    if strcmp(data.options.transcription,'multiple_shooting')


        solution.X=reshape(data.map.Vx*z,n,M)';
        solution.x0=solution.X(1,:)';
        usp=data.map.Vu*z;
        solution.U=reshape([usp;usp(end-m+1:end)],m,M)';
        solution.T=(solution.tf-t0)*[0;cumsum(data.tau)]+t0; 

    else

        solution.X=reshape(data.map.Vx*z,n,M)';
        solution.x0=solution.X(1,:)';
        usp=reshape(data.map.Vu*z,m,N)';
        solution.U=kron(usp,ones((M)/N,1));
        solution.T=(solution.tf-t0)*[0;cumsum(data.tau)*data.Nm/ns]+data.k0;

        solution.multipliers.lambdaNLP=solution.multipliers.lambda;

%         solution.multipliers.lambda_1toN=reshape(solution.multipliers.lambdaNLP(n+1:n*M)'*data.map.B,n,M)';
        solution.multipliers.lambda_1toM=reshape(solution.multipliers.lambdaNLP(1:n*M)',n,M)';
        
        if strcmp(data.options.transcription,'hermite')
            solution.X=solution.X(1:2:end,:);
            solution.U=solution.U(1:2:end,:);
            solution.T=solution.T(1:2:end);
            solution.multipliers.lambda_1toN=solution.multipliers.lambda_1toM(1:2:end,:);
        end
        

        
        if ng
            solution.multipliers.lambda_g=reshape(solution.multipliers.lambdaNLP(n*M+1:n*M+ng*M)',ng,M)';
        end
        if nrc
            solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP(n*M+ng*M+1:n*M+ng*M+nrc);
        end
        if nb
            solution.multipliers.lambda_b=solution.multipliers.lambdaNLP(n*M+M*ng+nrc+(~~nb):n*M+M*ng+nrc+nb);
        end
    end
end

%scaling back
if data.options.scaling
    solution.X=scale_variables_back( solution.X, data.data.Xscale, data.data.Xshift );
    solution.x0=scale_variables_back( solution.x0', data.data.Xscale, data.data.Xshift )';
    solution.U=scale_variables_back( solution.U, data.data.Uscale, data.data.Ushift );
    if isfield(data.data,'Pscale')
        solution.p=scale_variables_back( solution.p, data.data.Pscale, data.data.Pshift );
    end
    solution.tf=scale_variables_back( solution.tf, data.data.Tscale, data.data.Tshift );
    solution.T=scale_variables_back( solution.T, data.data.Tscale, data.data.Tshift );
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