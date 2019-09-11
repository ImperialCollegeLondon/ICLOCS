function [solution,XU0f]=resMinOptimization(required,z,data,phaseNo)
%resMinOptimization - Generate the cost, constraint and gradient
%information for the integrated residual minimization method
%
% Syntax:  solution=directCollocation(required,z,data)
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%    phaseNo - phase number
%
% Outputs:
%    solution - Data structure containing the solution
%    XU0f - data for linkage constraints
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

%%
global sol;

dataNLP=data.dataNLP;
[L,E,f,g,avrc,b]=deal(dataNLP.functions{:});
sol{phaseNo}.z=z;

% Get matrices
mp=dataNLP.map;

switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        [nt,np,n,m,ng,~,M,~,~,npd,~,~,nps,~,~,~,~]=deal(dataNLP.sizes{1:17});

        if np
            if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
                p=data.P;
                P=reshape(repmat(p,M,1),M,np);
            else
                p=z(n*(M+1)+M*m+1:n*(M+1)+M*m+np);
                P=reshape(repmat(p,M,1),M,np);
            end
        else
            p=[];
            P=spalloc(M,0,1);
        end

        if isfield(data,'free_time') && data.free_time
            Z=z;
            tf=z(end);
            t0=z(end-nt+1);
        else
            Z=[z;zeros(nt,1)];
            tf=data.tf_org;
            t0=data.t0_org;
        end

        X=reshape(mp.Vx*Z,M+1,n); %States including the end point
        U_Nm1=reshape(mp.Vu*Z,M,m); %input exluding the end point
        U=[U_Nm1;U_Nm1(end,:)];

        x0=X(1,:)';
        u0=U(1,:)';
        xf=X(end,:)';
        uf=U(end,:)';


        T=[data.tau;1]*(tf-t0)+t0;
        tau=[dataNLP.tau_inc;1];
        %for flexible mesh
        t_segment=(tf-t0)/2*dataNLP.tau_segment'+(tf+t0)/2; %Time at start/end of each segment
        t_segment_hes=[t_segment(1); t_segment(end)];

        t_segment_mat_m=data.dataNLP.t_segment_mat_m;
        t_segment_end=data.dataNLP.t_segment_mat_m*t_segment;
        data.dataNLP.t_segment_end=t_segment_end;

    otherwise % h Transcription Method

        [nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(dataNLP.sizes{1:13});

        if np
            if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
                p=data.P;
                P=reshape(repmat(p,M,1),M,np);
            else
                p=z(nt+1:nt+np);
                P=reshape(repmat(p,M,1),np,M)';
            end
        else
            p=[];
            P=spalloc(M,0,1);
        end

        if nt==1
            t0=dataNLP.t0;
            tf=z(1);
        elseif nt>=2
            t0=z(1);
            tf=z(2);
        else
            t0=dataNLP.t0;
            tf=dataNLP.tf; 
        end

        Z=z;

        X=reshape(mp.Vx*Z,n,M)';
        U=reshape(mp.Vu*Z,m,N)';

        x0=X(1,:)';
        u0=U(1,:)';
        xf=X(end,:)';
        uf=U(end,:)';

        tau=data.tau;
        T=tau*(tf-t0)+t0;


end

if strcmp(dataNLP.options.derivatives,'adigator')
    [ ResNorm_vec,Res_vec ] = getAdigator4ICLOCS_minres( X, U, T, P, data );
end

% Format reference inputs and states if applicable
if ~isempty(dataNLP.references.xr)
    Xr=dataNLP.references.xr;
else
    Xr=[];
end

if ~isempty(dataNLP.references.ur)
    Ur=dataNLP.references.ur;
else
    Ur=[];
end

% Prepare the data for linkage constraints
[XU0f] = prepareForLinkageConstraintInfo(x0,xf,u0,uf,data.dataNLP.data);


%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------
switch required

    case{'cost'}
        if data.mode==0
            switch dataNLP.options.discretization
                case{'globalLGR','hpLGR'} % p/hp Transcription Method
                    X_mesh=X(1:end-1,:);
                    U_mesh=U(1:end-1,:);
                    T_mesh=T(1:end-1,:);
                    Cost=data.dataNLP.map.w'*((tf-t0)*data.DT_seg_mesh./2.*L(X_mesh,Xr,U_mesh,Ur,P,T_mesh,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
                otherwise % h Transcription Method
                    Cost=data.dataNLP.map.w'*((tf-t0)*L(X,[],U,[],P,T,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
            end
            if strcmp(dataNLP.options.derivatives,'adigator')
                ResNorm_v=ResNorm_vec.f;
                sol{phaseNo}.cost= Cost+ResNorm_v;
            else
                [ResNorm,~] = costResidualMin_ModeMinRes( X,U,P,T,data);
                sol{phaseNo}.cost= Cost+sum(ResNorm,2);
            end
            
        elseif data.mode==1
            switch dataNLP.options.discretization
                case{'globalLGR','hpLGR'} % p/hp Transcription Method
                    X_mesh=X(1:end-1,:);
                    U_mesh=U(1:end-1,:);
                    T_mesh=T(1:end-1,:);
                    Cost=data.dataNLP.map.w'*((tf-t0)*data.DT_seg_mesh./2.*L(X_mesh,Xr,U_mesh,Ur,P,T_mesh,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
                otherwise % h Transcription Method
                    Cost=data.dataNLP.map.w'*((tf-t0)*L(X,[],U,[],P,T,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
            end
            sol{phaseNo}.cost= Cost;
        elseif data.mode==2
            if strcmp(dataNLP.options.derivatives,'adigator')
                ResNorm_v=ResNorm_vec.f;
                sol{phaseNo}.cost= ResNorm_v;
            else
                [ResNorm,~] = costResidualMin_ModeMinRes( X,U,P,T,data);
                sol{phaseNo}.cost= sum(ResNorm,2);
            end
        end
        solution=sol{phaseNo}.cost;


    case{'gradCost'}
        if data.mode==0
            switch dataNLP.options.discretization
                case{'globalLGR','hpLGR'} % p/hp Transcription Metho
                    X_mesh=X(1:end-1,:);
                    U_mesh=U(1:end-1,:);
                    T_mesh=T(1:end-1,:);
                    [gradCost,~]=gradientCostLGR(L,X_mesh,Xr,U_mesh,Ur,P,dataNLP.tau_inc,T_mesh,E,x0,xf,u0,uf,p,t0,tf,data.dataNLP);
                otherwise % h Transcription Method
                    [gradCost,~]=gradientCost(L,X,Xr,U,Ur,P,tau,E,x0,xf,u0,uf,p,t0,tf,data.dataNLP);
            end
            sol{phaseNo}.gradCost=gradCost+resMin_gradientCost_ModeMinRes(X,U,P,tau,t0,tf,data);
        elseif data.mode==1
            switch dataNLP.options.discretization
                case{'globalLGR','hpLGR'} % p/hp Transcription Metho
                    X_mesh=X(1:end-1,:);
                    U_mesh=U(1:end-1,:);
                    T_mesh=T(1:end-1,:);
                    [sol{phaseNo}.gradCost,sol{phaseNo}.JL]=gradientCostLGR(L,X_mesh,Xr,U_mesh,Ur,P,dataNLP.tau_inc,T_mesh,E,x0,xf,u0,uf,p,t0,tf,data.dataNLP);
                otherwise % h Transcription Method
                    [sol{phaseNo}.gradCost,sol{phaseNo}.JL]=gradientCost(L,X,Xr,U,Ur,P,tau,E,x0,xf,u0,uf,p,t0,tf,data.dataNLP);
            end
        elseif data.mode==2
            if strcmp(dataNLP.options.derivatives,'adigator')
                sol{phaseNo}.gradCost=zeros(1,ResNorm_vec.dY_size);
                sol{phaseNo}.gradCost(1,ResNorm_vec.dY_location)=ResNorm_vec.dY;
            else
                sol{phaseNo}.gradCost=resMin_gradientCost_ModeMinRes(X,U,P,tau,t0,tf,data);
            end
            
        end
        solution=sol{phaseNo}.gradCost; 


    case{'const'}
        
        if data.mode==1 
            if strcmp(dataNLP.options.derivatives,'adigator')
                const = constResidualMin_ModeMinCost_Adigator( X, U, P, T, data, 1);
                const= [const;Res_vec.f];
            else
                [ const,~ ] = constResidualMin_ModeMinCost( X, U, P, T, data, 1);
            end
        else
            if strcmp(dataNLP.options.derivatives,'adigator')
                Res=Res_vec.f;
                sol{phaseNo}.constRes = Res;
                const = [constResidualMin_ModeMinRes( X, U, P, T, data);Res];
            else
                [~,Res]=costResidualMin_ModeMinRes( X,U,P,T,data);
                Res=sum(Res,2);
                sol{phaseNo}.constRes = Res;
                const = [constResidualMin_ModeMinRes( X, U, P, T, data);Res];
            end
        end
        sol{phaseNo}.const = const;
        solution=sol{phaseNo}.const;
    case{'jacConst'}
        if data.mode==1
            if strcmp(dataNLP.options.derivatives,'adigator')
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        jac=resMin_jacobianFD_LGR_ModeMinCost_Adigator(g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,data,data.free_time,Res_vec);
                    otherwise
                        jac=resMin_jacobianFD_ModeMinCost_Adigator(g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,data,Res_vec);
                end
            else
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        jac=resMin_jacobianFD_LGR_ModeMinCost(g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,data,data.free_time);
                    otherwise
                        jac=resMin_jacobianFD_ModeMinCost(g,avrc,X,U,P,tau,b,x0,xf,u0,uf,p,t0,tf,data);
                end
            end

        else
            if strcmp(dataNLP.options.derivatives,'adigator')
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        jac=resMin_jacobianFD_LGR_ModeMinRes_Adigator(L,E,g,avrc,X,Xr,U,Ur,P,tau,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_mat_m,data,data.free_time,Res_vec);
                    otherwise
                        jac=resMin_jacobianFD_ModeMinRes_Adigator(L,E,g,avrc,X,Xr,U,Ur,P,tau,b,x0,xf,u0,uf,p,t0,tf,data,Res_vec);
                end
            else
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        jac=resMin_jacobianFD_LGR_ModeMinRes(L,E,g,avrc,X,Xr,U,Ur,P,tau,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_mat_m,data,data.free_time);
                    otherwise
                        jac=resMin_jacobianFD_ModeMinRes(L,E,g,avrc,X,Xr,U,Ur,P,tau,b,x0,xf,u0,uf,p,t0,tf,data);
                end
            end
            

        end
        sol{phaseNo}.jacConst=sparse(jac);
        solution=sol{phaseNo}.jacConst;

    case{'hessian'}
        if data.mode==0 %Weighted
             switch dataNLP.options.discretization
                case{'globalLGR','hpLGR'} % p/hp Transcription Method
                    Hes=resMin_hessianCD_LGR_ModeMinCost(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_hes,data);
                otherwise
                    Hes=resMin_hessianCD_ModeMinWeightedCostRes(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data);
            end
        elseif data.mode==1 %Minimizing Cost
            if strcmp(dataNLP.options.derivatives,'adigator')
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        Hes=resMin_hessianCD_LGR_ModeMinCost_Adigator(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_hes,data,Res_vec);
                    otherwise
                        Hes=resMin_hessianCD_ModeMinCost_Adigator(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data,Res_vec);
                end
            else
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        Hes=resMin_hessianCD_LGR_ModeMinCost(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_hes,data);
                    otherwise
                        Hes=resMin_hessianCD_ModeMinCost(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data);
                end
            end

        elseif data.mode==2 %Minimizing Residual
            if strcmp(dataNLP.options.derivatives,'adigator')
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        Hes=resMin_hessianCD_LGR_ModeMinRes_Adigator(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_hes,data,ResNorm_vec,Res_vec);
                    otherwise
                        Hes=resMin_hessianCD_ModeMinRes_Adigator(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data,ResNorm_vec,Res_vec);
                end
            else
                switch dataNLP.options.discretization
                    case{'globalLGR','hpLGR'} % p/hp Transcription Method
                        Hes=resMin_hessianCD_LGR_ModeMinRes(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment_hes,data);
                    otherwise
                        Hes=resMin_hessianCD_ModeMinRes(L,g,X,U,P,tau,E,b,x0,xf,u0,uf,p,t0,tf,data);
                end
            end

        end
        sol{phaseNo}.hessian=sparse(tril(Hes));
        solution=sol{phaseNo}.hessian;
end

