function [ varargout ] = solInterpolationMinResidual( problem,solution,t_list,options,dataNLP )
%solInterpolationMinResidual - solution representation with integrated
%residual minimization
%
% Syntax:   [X_quad,dX_quad,U_quad,tau_quad,LGR_points_quad,npd_quad,t0,tf,p,scaled_solution,D_mat] = solInterpolationMinResidual( problem,solution,t_list,options,dataNLP) - hpLGR method
%           [Xp,dXp,Up,t0,tf,p,minres_solution,scaled_solution] = solInterpolationMinResidual( problem,solution,t_list,options,dataNLP) - Hermite-Simpson method
% 
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------


[nt,np,n,m,ng,~,M]=deal(dataNLP.sizes{1:7});

auxdata = transcribeResMin( t_list,options,dataNLP );

        
auxdata.cost_org=solution.cost;
auxdata.tf_org=solution.tf;

switch dataNLP.options.discretization
    
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        idx_endpoint=dataNLP.tau_local_seg == -1;
        if dataNLP.options.scaling
            auxdata.dataNLP.data.u_endpoint=[solution.scaledVariables.U(idx_endpoint,:);solution.scaledVariables.U(end,:)];
            auxdata.dataNLP.data.x_endpoint=[solution.scaledVariables.X(idx_endpoint,:);solution.scaledVariables.X(end,:)];  
        else
            auxdata.dataNLP.data.u_endpoint=[solution.U(idx_endpoint,:);solution.U(end,:)];
            auxdata.dataNLP.data.x_endpoint=[solution.X(idx_endpoint,:);solution.X(end,:)];  
        end


        if ~strcmp(options.derivatives,'adigator')
            auxdata.free_time=~(dataNLP.infoNLP.zu(end)==dataNLP.infoNLP.zl(end));
        else
            auxdata.free_time=1;
        end
        
        if auxdata.free_time
            options.lb = dataNLP.infoNLP.zl;  % Lower bound on the variables.
            options.ub = dataNLP.infoNLP.zu;  % Upper bound on the variables.
            z0     = solution.z;
            options.zl     = dataNLP.infoNLP.zl;
            options.zu     = dataNLP.infoNLP.zu;
        else
            options.lb = dataNLP.infoNLP.zl(1:(end-nt));  % Lower bound on the variables.
            options.ub = dataNLP.infoNLP.zu(1:(end-nt));  % Upper bound on the variables.
            z0     = solution.z(1:(end-nt));
            options.zl     = dataNLP.infoNLP.zl(1:(end-nt));
            options.zu     = dataNLP.infoNLP.zu(1:(end-nt));
        end

        auxdata.idx_endpoint=idx_endpoint;
        
        auxdata.t0_org=solution.T(1);
        
        if dataNLP.options.resminRep.collmatch
            n_const=2*(n);
        else
            n_const=0;
        end
       
    otherwise % h Transcription Method
        if dataNLP.options.scaling
            auxdata.dataNLP.data.u_endpoint=solution.scaledVariables.U;
            auxdata.dataNLP.data.x_endpoint=solution.scaledVariables.X;
        else
            auxdata.dataNLP.data.u_endpoint=solution.U;
            auxdata.dataNLP.data.x_endpoint=solution.X;
        end
        
        if nt
            options.lb = dataNLP.infoNLP.zl;  % Lower bound on the variables.
            options.ub = dataNLP.infoNLP.zu;  % Upper bound on the variables.
            z0     = solution.z;
            options.zl     = dataNLP.infoNLP.zl;
            options.zu     = dataNLP.infoNLP.zu;
        else
            options.lb = dataNLP.infoNLP.zl((1+nt):end);  % Lower bound on the variables.
            options.ub = dataNLP.infoNLP.zu((1+nt):end);  % Upper bound on the variables.
            z0     = solution.z((1+nt):end);
            options.zl     = dataNLP.infoNLP.zl((1+nt):end);
            options.zu     = dataNLP.infoNLP.zu((1+nt):end);
        end
        auxdata.t_mesh=[0;dataNLP.tau_inc]/2;
        auxdata.t0_org=solution.t0;
        
        if dataNLP.options.resminRep.collmatch
            n_const=2*(n);
        else
            n_const=0;
        end
end

if dataNLP.options.scaling
    auxdata.P=solution.scaledVariables.p;
else
    auxdata.P=solution.p;
end

auxdata.dataNLP.data.discErrorConstScaling=ones(1,n);

% formulate constraints
switch dataNLP.options.discretization
    
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        if isfield(solution,'scaledVariables')
            [~,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,[solution.scaledVariables.coll.U;solution.scaledVariables.coll.U(end,:)],reshape(repmat(solution.scaledVariables.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],auxdata);
        else
            [~,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,[solution.coll.U;solution.coll.U(end,:)],reshape(repmat(solution.coll.p,M+1,1),M+1,np),[solution.coll.T;solution.tf],auxdata);
        end
        ResColl=sum(ResColl,2);
        if  dataNLP.options.resminRep.collmatch
            options.cl = [dataNLP.glAll;dataNLP.map.rcl(:);-eps*ones(n_const,1);-inf;zeros(length(ResColl),1)];   % Lower bounds on the constraint functions.
            options.cu = [dataNLP.guAll;dataNLP.map.rcu(:);eps*ones(n_const,1);auxdata.cost_org;ones(length(ResColl),1)];  % Upper bounds on the constraint functions.
            c_Tol=[dataNLP.g_tolAll;zeros(size(dataNLP.map.rcu(:)));zeros(n_const,1);0;zeros(size(ResColl))];
        else
            options.cl = [dataNLP.glAll;dataNLP.map.rcl(:);dataNLP.map.bl(:);-inf;zeros(length(ResColl),1)];   % Lower bounds on the constraint functions.
            options.cu = [dataNLP.guAll;dataNLP.map.rcu(:);dataNLP.map.bu(:);auxdata.cost_org;ones(length(ResColl),1)];  % Upper bounds on the constraint functions.
            c_Tol=[dataNLP.g_tolAll;zeros(size(dataNLP.map.rcu(:)));dataNLP.map.bTol(:);0;zeros(size(ResColl))];
        end

    otherwise % h Transcription Method
        if isfield(solution,'scaledVariables')
            [~,ResColl]=costResidualMin_ModeMinRes( solution.scaledVariables.coll.X,solution.scaledVariables.coll.U,reshape(repmat(solution.scaledVariables.coll.p,M,1),M,np),solution.coll.T,auxdata);
        else
            [~,ResColl]=costResidualMin_ModeMinRes( solution.coll.X,solution.coll.U,reshape(repmat(solution.coll.p,M,1),M,np),solution.coll.T,auxdata);
        end
        ResColl=sum(ResColl,2);
        options.ipopt.constr_viol_tol=min(ResColl);
        if  dataNLP.options.resminRep.collmatch
            options.cl = [dataNLP.glAll;dataNLP.map.rcl(:);-eps*ones(n_const,1);-inf;zeros(length(ResColl),1)];   % Lower bounds on the constraint functions.
            options.cu = [dataNLP.guAll;dataNLP.map.rcu(:);eps*ones(n_const,1);auxdata.cost_org;ones(length(ResColl),1)];  % Upper bounds on the constraint functions.
            c_Tol = [dataNLP.g_tolAll;zeros(size(dataNLP.map.rcu(:)));zeros(n_const,1);0;zeros(size(ResColl))];  % Upper bounds on the constraint functions.
        else
            options.cl = [dataNLP.glAll;dataNLP.map.rcl(:);dataNLP.map.bl(:);-eps*ones(n_const,1);-inf;zeros(length(ResColl),1)];   % Lower bounds on the constraint functions.
            options.cu = [dataNLP.guAll;dataNLP.map.rcu(:);dataNLP.map.bu(:);eps*ones(n_const,1);auxdata.cost_org;ones(length(ResColl),1)];  % Upper bounds on the constraint functions.
            c_Tol = [dataNLP.g_tolAll;zeros(size(dataNLP.map.rcu(:)));dataNLP.map.bTol(:);0;zeros(size(ResColl))];  % Upper bounds on the constraint functions.
        end
end
auxdata.infoNLP.cu=options.cu;
auxdata.infoNLP.cl=options.cl;
auxdata.infoNLP.zu=options.zu;
auxdata.infoNLP.zl=options.zl;
auxdata.infoNLP.c_Tol=c_Tol;
auxdata.dataNLP.data.discErrorConstScaling=1./ResColl';
ResConstScaleMat=repmat(auxdata.dataNLP.data.discErrorConstScaling, auxdata.nps, 1 );
auxdata.ResConstScaleMat=diag(ResConstScaleMat(:));

        
discErrorTol_Full=problem.states.xErrorTol_integral.^2;
if any(discErrorTol_Full<eps)
    error('Integral of the residual errors allowed must be strictly larger than machine precision');
else
    discErrorTol_Full=discErrorTol_Full';
end
auxdata.dataNLP.data.discErrorTol_Full=discErrorTol_Full;
auxdata.dataNLP.data.discErrorTol_FullScaling=discErrorTol_Full.*auxdata.dataNLP.data.discErrorConstScaling';
 
n_const=length(options.cl);
auxdata.n_const=n_const;
  
% Initialize the dual point.
options.lambda = zeros(n_const,1);
  
% Set the IPOPT options.
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.print_level = 5;
options.ipopt.tol         = dataNLP.options.ipopt.tol;
%   options.ipopt.hessian_approximation='limited-memory';
options.ipopt.hessian_approximation='exact';
  
% The callback functions.
auxdata.options.transcription='direct_collocation';
auxdata.options.discretization='resMinInterpolationForSolution';

auxdata.mode=2;
if strcmp(options.derivatives,'adigator')
    genAdigator4ICLOCS_resmin( options, auxdata, n, m, np, nt, M );
end

    
options.auxdata=auxdata;
data.funcs.objective         = @costFunction;
data.funcs.gradient          = @costGradient;
data.funcs.constraints       = @constraintFunction;
data.funcs.jacobian          = @constraintJacobian;
data.funcs.jacobianstructure = @(data) sparse(ones(n_const,length(z0)));
data.funcs.hessian           = @computeHessian;
data.funcs.hessianstructure  = @(data) sparse(tril(ones(length(z0),length(z0))));

if options.resminEarlyStop
    data.funcs.iterfunc=@callback_minres;
end

% Run IPOPT.
[z,~] = ipopt_auxdata(z0,data.funcs,options);

switch dataNLP.options.discretization
    
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        [nt,np,n,m,ng,~,M,~,~,~,~,~,~,~,~,~,~]=deal(dataNLP.sizes{1:17});
        if auxdata.free_time
            tf=z(end);
            t0=z(end-nt+1);
        else
            z=[z;zeros(nt,1)];
            t0=auxdata.t0_org;
            tf=auxdata.tf_org;
        end
        
        X_Np1=reshape(dataNLP.map.Vx*z,M+1,n); %States including the end point
        X=X_Np1(1:M,:); %States exlduing the end point
        U=reshape(dataNLP.map.Vu*z,M,m); %input exluding the end point
        U_Np1=[U;U(end,:)];
        
        if dataNLP.options.resminRep.collmatch
            p=solution.p;
        else
            p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np);
        end
        
        X_quad=(auxdata.InterpH*X_Np1)./auxdata.sumInterpH;
        U_quad=(auxdata.InterpH*U_Np1)./auxdata.sumInterpH;
        X_quad(auxdata.interp_fixi,:)=X_Np1(auxdata.interp_fixj,:);
        U_quad(auxdata.interp_fixi,:)=U_Np1(auxdata.interp_fixj,:);

        DT=(auxdata.tf_list-auxdata.t0_list)/2;
        DT_seg_quad=repelem(DT',auxdata.npd_quad+1,1);
        dX_quad=auxdata.D_mat*X_quad./DT_seg_quad/(tf-t0);

        LGR_points_quad=[auxdata.LGR.points{1,auxdata.idx_quad};1];
        
        if dataNLP.options.scaling
            scaled_solution.X=X;
            scaled_solution.U=U;
            X_quad=scale_variables_back( X_quad, dataNLP.data.Xscale, dataNLP.data.Xshift );
            U_quad=scale_variables_back( U_quad, dataNLP.data.Uscale, dataNLP.data.Ushift );
            dX_quad=scale_variables_back( dX_quad, dataNLP.data.Xscale, 0 );
            if isfield(dataNLP.data,'Pscale')
                scaled_solution.p=p;
                p=scale_variables_back( p', dataNLP.data.Pscale, dataNLP.data.Pshift )';
            end
        else
            scaled_solution=[];
        end
        
        varargout{1}=X_quad;
        varargout{2}=dX_quad;
        varargout{3}=U_quad;
        varargout{4}=auxdata.tau_quad;
        varargout{5}=LGR_points_quad;
        varargout{6}=auxdata.npd_quad;
        varargout{7}=t0;
        varargout{8}=tf;
        varargout{9}=p;
        varargout{10}=scaled_solution;
        varargout{11}=auxdata.D_mat;
    case{'hermite'} % h Transcription Method
        [nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(dataNLP.sizes{1:13});
        if nt==1
            t0=auxdata.t0_org;
            tf=z(1);
        elseif nt==2
            t0=z(1);
            tf=z(2);
        else
            z=[zeros(nt,1);z];
            t0=auxdata.t0_org;
            tf=auxdata.tf_org;
        end
        X=reshape(dataNLP.map.Vx*z,n,M)';
        U=reshape(dataNLP.map.Vu*z,m,N)';
        
        if dataNLP.options.resminRep.collmatch
           p=solution.p;
        else
           p=z(nt+1:nt+np);
        end

        f=dataNLP.functions{3};
        F_k=f(X(1:2:end,:),U(1:2:end,:),repmat(p,length(auxdata.tau),1),auxdata.tau(1:2:end)*(tf-t0),dataNLP.data);
        F_kph=auxdata.DxHS_hf*X/(tf-t0)-F_k(1:end-1,:)/2;
        F_kp1=auxdata.DxHS_p1*X/(tf-t0)+F_k(1:end-1,:);
        F=[F_k(1:end-1,:) F_kph F_kp1]';
        F=reshape(F(:),n,3*auxdata.nps)';

        T=auxdata.tau*(tf-t0)+t0;

        %scaling back
        if dataNLP.options.scaling
            scaled_solution.X=X;
            scaled_solution.U=U;
            X=scale_variables_back( X, dataNLP.data.Xscale, dataNLP.data.Xshift );
            U=scale_variables_back( U, dataNLP.data.Uscale, dataNLP.data.Ushift );
            F=scale_variables_back( F, dataNLP.data.Xscale, 0 );
            if isfield(dataNLP.data,'Pscale')
                scaled_solution.p=p;
                p=scale_variables_back( p', dataNLP.data.Pscale, dataNLP.data.Pshift )';
            end
        else
            scaled_solution=[];
        end

        minres_solution.X=X;
        minres_solution.U=U;
        minres_solution.p=p;
        minres_solution.T=T;
        
        for i=1:n % Cubic Hermite interpolation
            [Xp{i}, dXp{i}]=HSInterpolation(minres_solution.T,minres_solution.X(:,i),F(:,i));
        end
        for i=1:m % Piecewise quadratic polynomials
            [Up{i}]=HSInterpolationU(minres_solution.T,minres_solution.U(:,i));
        end
        varargout{1}=Xp;
        varargout{2}=dXp;
        varargout{3}=Up;
        varargout{4}=t0;
        varargout{5}=tf;
        varargout{6}=p;
        varargout{7}=minres_solution;
        varargout{8}=scaled_solution;
        
    otherwise
        [nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(dataNLP.sizes{:});
        z=[zeros(np+nt,1);z];
        X=reshape(dataNLP.map.Vx*z,n,M)';
        U=reshape(dataNLP.map.Vu*z,m,N)';
        for i=1:n % Cubic Hermite interpolation
            [Xp{i}, dXp{i}]=Hsplines(T,X(:,i),F(:,i));
        end
        for i=1:m % Piecewise quadratic polynomials
            [Up{i}]=QuadsplinesU(T,U(:,i));
        end
end




end

