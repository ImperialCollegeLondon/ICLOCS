function [solution,status,data]=solveSingleNLP_ResidualMin(NLP,data)

%solveSingleNLP_DirectCollocation - Solve a single nonlinear program and
%return the solution (integrated residual minimization method)
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

data.infoNLP=NLP;

if strcmp(data.dataNLP.options.discretization,'globalLGR') || strcmp(data.dataNLP.options.discretization,'hpLGR')
    [nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.dataNLP.sizes{1:17});
    if data.dataNLP.options.reorderLGR
        NLP.z0(1:n)=data.dataNLP.x0t;
    else
        NLP.z0(1:M+1:M*n)=data.dataNLP.x0t;
    end
else
    [nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.dataNLP.sizes{1:13});
    NLP.z0(nt+np+1:n+nt+np)=data.dataNLP.x0t;
end
nrc=nrcl+nrcu+nrce;

  

switch(data.dataNLP.options.NLPsolver)

    case{'ipopt'} 
       if ~isempty(data.dataNLP.multipliers)
            opt.lambda=data.dataNLP.multipliers.lambda;              % Solve the NLP using IPOPT
            data.dataNLP.options.ipopt.warm_start_init_point ='yes'; 
       end
       opt.ipopt=data.dataNLP.options.ipopt;
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
       
       if  ((strcmp(data.dataNLP.options.discretization,'globalLGR')) || (strcmp(data.dataNLP.options.discretization,'hpLGR'))) && data.dataNLP.options.reorderLGR
           solution.z_org=z;
           z=z(data.dataNLP.reorder.z_idx_back);
           info.zl=info.zl(data.dataNLP.reorder.z_idx_back);
           info.zu=info.zu(data.dataNLP.reorder.z_idx_back);
           info.lambda=info.lambda(data.dataNLP.reorder.vert_idx_back);
       end
       
       solution.multipliers.zl=info.zl;
       solution.multipliers.zu=info.zu;
       solution.multipliers.lambda=info.lambda;
       status=info;
       
    
    case{'fmincon'}                               % Solve the NLP using FMINCON
        error('use of fmincon not yet supported for direct integrated residual minimization method, please use IPOPT')
        prompt = 'WARNING: Use of fmincon is NOT ADVISABLE because the solution speed can be significantly slower with a higher chance of failure. \n We would recommend to use fmincon only for sanity checks. \n For solving, please consider use IPOPT (see http://www.ee.ic.ac.uk/ICLOCS/Downloads.html). \n Press any key to continue...\n ';
        input(prompt)
        tA=tic;
        [z,cost,status,output,multipliers]=fmincon(@(z)fminCost(z,data),NLP.z0,...
            [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.fmincon);
        tB=toc(tA);
        if  ((strcmp(data.dataNLP.options.discretization,'globalLGR')) || (strcmp(data.dataNLP.options.discretization,'hpLGR'))) && data.dataNLP.options.reorderLGR
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
        if isfield(data.dataNLP.options,'NOMAD')
           opts = nomadset('display_degree',2,'bb_output_type',obj_type,'initial_mesh_size','50000','h_min',1e-03,'f_target',data.options.NOMAD.f_target);
        else
%            opts = nomadset('display_degree',2,'bb_output_type',obj_type,'initial_mesh_size','50000','h_min',1e-03,'h_norm','L2');
           opts = nomadset('display_degree',2,'bb_output_type',obj_type);

        end
        xtype    = cell(length(NLP.z0),1);
        xtype(:) = {'C'};
        xtype=char(xtype);
        [z,cost,status,output.iter,output.nfval] = nomad(@(z)MADSCost(z,data),NLP.z0,NLP.zl,NLP.zu,@(z)MADSConst(z,data,NLP),[NLP.cu;-NLP.cl],xtype,opts);
%         [z,cost,status,output.iter,output.nfval] = nomad(@(z)MADSCost(z,data),NLP.z0,NLP.zl,NLP.zu,[],[],xtype,opts);

%         [z,cost,status,output]=patternsearch(@(z)fminCost(z,data),NLP.z0,...
%             [],[],[],[],NLP.zl,NLP.zu,@(z)fminConst(z,data,NLP),data.options.MADS);
        tB=toc(tA);
        if  ((strcmp(data.dataNLP.options.discretization,'globalLGR')) || (strcmp(data.dataNLP.options.discretization,'hpLGR'))) && data.dataNLP.options.reorderLGR
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
solution.computation_time=tB;

% solution.z=NLP.z0;
solution.z=z;


if strcmp(data.dataNLP.options.discretization,'globalLGR') || strcmp(data.dataNLP.options.discretization,'hpLGR')
    if nt;solution.tf=z(end); else solution.tf=data.dataNLP.tf; end
    if np;solution.coll.p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np); else solution.coll.p=[]; end
    t0=z(end-nt+1);
    solution.t0=t0;

    solution.coll.tf=solution.tf;
    solution.coll.X=reshape(data.dataNLP.map.Vx*z,M+1,n);
    solution.coll.x0=solution.coll.X(1,:)';
    solution.coll.U=reshape(data.dataNLP.map.Vu*z,M,m);
    
    if data.dataNLP.options.adaptseg==1
        t_segment=z(end-nt+1:end)';

        tau_segment = ((1-(-1))*t_segment-(t_segment(1)+t_segment(end))) /(t_segment(end)-t_segment(1));

        tau_inc=[];

        tau_seg=cell(1,length(tau_segment)-1);
        for i=1:length(tau_segment)-1
            tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [data.dataNLP.map.LGR.points{npduidx(i)}; 1]) / (1-(-1));
            tau_inc=[tau_inc; tau_seg{i}(1:end-1)];
        end

        tau=diff(tau_inc);
        data.dataNLP.tau_inc=tau_inc;
        data.dataNLP.tau_segment=tau_segment;
        data.dataNLP.tau=tau;
    end
    
    solution.coll.T=(solution.tf-t0)/2*data.dataNLP.tau_inc+(solution.tf+t0)/2; 
    
    if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 
        solution.multipliers.lambdaNLP=solution.multipliers.lambda;
%         solution.multipliers.lambda_1toN=reshape(solution.multipliers.lambda(1:M*n),M,n);
%         D_N=diag(1./data.dataNLP.map.w);
%         solution.multipliers.lambda=D_N*solution.multipliers.lambda_1toN;
%         D_Np1=data.dataNLP.map.LGR.diff_matrix{end}(:,end);
%         solution.multipliers.lambda(end+1,:)=D_Np1'*solution.multipliers.lambda_1toN(end-length(data.dataNLP.map.LGR.points{end})+1:end,:);
        if ng
            solution.multipliers.lambda_g=zeros(ng,M);
            solution.multipliers.lambda_g(logical(data.dataNLP.gActiveIdx'))=solution.multipliers.lambdaNLP(1:ngActive);
            solution.multipliers.lambda_g=solution.multipliers.lambda_g';
        end
        if nrc
            solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP((ngActive+1):ngActive+nrc);
        end
        if nb
            solution.multipliers.lambda_b=solution.multipliers.lambdaNLP((ngActive+nrc+1):(ngActive+nrc+nb));
        end
        solution.multipliers.lambda_resconst=solution.multipliers.lambdaNLP(end-n+1:end);
    end
    
        solution.X=solution.coll.X;
        solution.U=solution.coll.U;
        solution.T=solution.coll.T;
        solution.p=solution.coll.p;
        solution.x0=solution.coll.x0;

else
    
    if nt==1
        solution.tf=z(1);
        t0=data.dataNLP.t0;
    elseif nt>=2
        solution.tf=z(2);
        t0=z(1);
    else
        solution.tf=data.dataNLP.tf; 
        t0=data.dataNLP.t0;
    end
    solution.t0=t0;
    
    if np
        solution.coll.p=z(nt+1:nt+np);
    else
        solution.coll.p=[]; 
    end
    

    if strcmp(data.dataNLP.options.transcription,'multiple_shooting')


        solution.X=reshape(data.dataNLP.map.Vx*z,n,M)';
        solution.x0=solution.X(1,:)';
        usp=data.dataNLP.map.Vu*z;
        solution.U=reshape([usp;usp(end-m+1:end)],m,M)';
        solution.T=(solution.tf-t0)*[0;cumsum(data.dataNLP.tau)]+t0; 

    else

        solution.coll.X=reshape(data.dataNLP.map.Vx*z,n,M)';
        solution.coll.x0=solution.coll.X(1,:)';
        usp=reshape(data.dataNLP.map.Vu*z,m,N)';
        solution.coll.U=kron(usp,ones((M)/N,1));
        solution.coll.T=(solution.tf-t0)*[0;cumsum(data.dataNLP.tau)*data.dataNLP.Nm/ns]+t0;

        if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 

            solution.multipliers.lambdaNLP=solution.multipliers.lambda;

    %         solution.multipliers.lambda_1toN=reshape(solution.multipliers.lambdaNLP(n+1:n*M)'*data.dataNLP.map.B,n,M)';
%             solution.multipliers.lambda_1toM=reshape(solution.multipliers.lambdaNLP(1:n*M)',n,M)';
        end
        
        if strcmp(data.dataNLP.options.discretization,'hermite')
            solution.X=solution.coll.X(1:2:end,:);
            solution.U=solution.coll.U(1:2:end,:);
            solution.T=solution.coll.T(1:2:end);
            if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 
%                 solution.multipliers.lambda_1toN=solution.multipliers.lambda_1toM(1:2:end,:);
            end
        else
            solution.X=solution.coll.X;
            solution.U=solution.coll.U;
            solution.T=solution.coll.T;
        end
        solution.p=solution.coll.p;
        solution.x0=solution.coll.x0;

        if ~strcmp(data.dataNLP.options.NLPsolver,'NOMAD') 
            if ng
                solution.multipliers.lambda_g=zeros(ng,M);
                solution.multipliers.lambda_g(logical(data.dataNLP.gActiveIdx'))=solution.multipliers.lambdaNLP(1:ngActive);
                solution.multipliers.lambda_g=solution.multipliers.lambda_g';
            end
            if nrc
                solution.multipliers.lambda_rc=solution.multipliers.lambdaNLP(ngActive+1:ngActive+nrc);
            end
            if nb
                solution.multipliers.lambda_b=solution.multipliers.lambdaNLP(ngActive+nrc+1:ngActive+nrc+nb);
            end
            solution.multipliers.lambda_resconst=solution.multipliers.lambdaNLP(end-n+1:end);
        end
    end
end

%scaling back
if data.dataNLP.options.scaling
    solution.scaledVariables.X=solution.X;
    solution.scaledVariables.x0=solution.x0;
    solution.scaledVariables.U=solution.U;
    solution.scaledVariables.coll.X=solution.coll.X;
    solution.scaledVariables.coll.x0=solution.coll.x0;
    solution.scaledVariables.coll.U=solution.coll.U;
    solution.X=scale_variables_back( solution.X, data.dataNLP.data.Xscale, data.dataNLP.data.Xshift );
    solution.x0=scale_variables_back( solution.x0', data.dataNLP.data.Xscale, data.dataNLP.data.Xshift )';
    solution.U=scale_variables_back( solution.U, data.dataNLP.data.Uscale, data.dataNLP.data.Ushift );
    solution.coll.X=scale_variables_back( solution.coll.X, data.dataNLP.data.Xscale, data.dataNLP.data.Xshift );
    solution.coll.x0=scale_variables_back( solution.coll.x0', data.dataNLP.data.Xscale, data.dataNLP.data.Xshift )';
    solution.coll.U=scale_variables_back( solution.coll.U, data.dataNLP.data.Uscale, data.dataNLP.data.Ushift );
    if isfield(data.dataNLP.data,'Pscale')
        solution.scaledVariables.p=solution.p;
        solution.scaledVariables.coll.p=solution.coll.p;
        solution.p=scale_variables_back( solution.p', data.dataNLP.data.Pscale, data.dataNLP.data.Pshift )';
    else
        solution.scaledVariables.p=[];
        solution.scaledVariables.coll.p=[];
    end
%     solution.tf=scale_variables_back( solution.tf, data.dataNLP.data.Tscale, data.dataNLP.data.Tshift );
%     solution.T=scale_variables_back( solution.T, data.dataNLP.data.Tscale, data.dataNLP.data.Tshift );
end


% if isfield(data.data,'singular_arc_lift')
%     solution.U(:,data.data.singular_arc_lift)=(solution.U(:,data.data.singular_arc_lift)+1).^2+data.data.singular_arc_lift_shift(data.data.singular_arc_lift);
% end
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