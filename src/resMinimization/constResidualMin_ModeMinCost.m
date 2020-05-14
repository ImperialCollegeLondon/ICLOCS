function [ const,Res ] = constResidualMin_ModeMinCost( X, U, P, T, data, mode)
%constResidualMin_ModeMinCost - constraint function computation for
%integrated residual minimization (alternating method: cost
%minimization)
%
% Syntax:   [ const,Res ] = constResidualMin_ModeMinCost( X, U, P, T, data, mode)
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

X_s=X;
U_s=U;
P_s=P;

if isfield(data.dataNLP.data,'Xscale')
    X=scale_variables_back( X, data.dataNLP.data.Xscale, data.dataNLP.data.Xshift );
    U=scale_variables_back( U, data.dataNLP.data.Uscale, data.dataNLP.data.Ushift );
    if isfield(data.dataNLP.data,'Pscale')
        P=scale_variables_back( P, data.dataNLP.data.Pscale, data.dataNLP.data.Pshift );
    end
end

dataNLP=data.dataNLP;

t0=T(1);tf=T(end);
x0=X(1,:);
u0=U(1,:)';
xf=X(end,:)';
uf=U(end,:)';

[~,~,f,g,avrc,b]=deal(data.dataNLP.functions_unscaled{:});

switch dataNLP.options.discretization
    
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        [~,~,~,~,ng,~,M,~,~,~,~,~,~,~,~,~,~,ng_eq,ng_neq]=deal(dataNLP.sizes{1:19});
        
        X_quad=(data.InterpH*X)./data.sumInterpH;
        U_quad=(data.InterpH*U)./data.sumInterpH;
        X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
        U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
        dX_quad=data.D_mat*X_quad./(data.DT_seg/2)/(tf-t0);
        P_quad=repelem(P(1,:),length(data.tau_quad),1);
        
        Fp=f(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);
        if ng_eq
            Gp_eq=g_eq(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);
            Gp_eq=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Gp_eq.^2));
        else
            Gp_eq=[];
        end
        Res=dX_quad-Fp;
        Res=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Res.^2));
        Res=scale_variables( [Res;Gp_eq]', data.dataNLP.data.discErrorConstScaling, 0 )';
        
        if mode==1
            X_mesh=X(1:end-1,:);
            U_mesh=U(1:end-1,:);
            T_mesh=T(1:end-1,:);
            X_mesh_Np1=X;
            
%             if ng_neq
%                 Gc=g(X_mesh,U_mesh,P,T_mesh,data.dataNLP.data);
%                 Gc_neq=Gc(:,ng_eq+1:end);
%             else
%                 Gc_neq=[];
%             end
%             g_vect=reshape(Gc_neq',M*ng,1);
            g_vect=reshape(g(X_mesh,U_mesh,P,T_mesh,data.dataNLP.data)',M*ng,1);
            if isfield(data.dataNLP.data,'Xscale')
                cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
            else
                cr=avrc(X_mesh_Np1,U_mesh,P,T,data.dataNLP)';
            end
            const=[  g_vect(data.dataNLP.gAllidx);
                     cr;
                     b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{});
                     sum(Res,2)];

        else
            const=[];
        end
                
    otherwise
        [~,~,n,~,ng,~,M,~,~,~,~,~,~,~,ng_eq,ng_neq]=deal(dataNLP.sizes{1:16});
        
        F_k=f(X(1:2:end,:),U(1:2:end,:),P(1:2:end,:),data.tau(1:2:end)*(tf-t0),data.dataNLP.data);
        F_kph=data.DxHS_hf*X/(tf-t0)-F_k(1:end-1,:)/2;
        F_kp1=data.DxHS_p1*X/(tf-t0)+F_k(1:end-1,:);
        F=[F_k(1:end-1,:) F_kph F_kp1]';
        F=reshape(F(:),n,3*data.nps)';
        
        X_quad=repelem(X(1:2:end-1,:),data.npd_quad+1,1)+(tf-t0)*data.AxHS*F;
        U_quad=data.AuHS*U;
        dX_quad=data.AfHS*F;
        
        P_quad=repelem(P(1,:),length(data.tau_quad),1);
        
        Fp=f(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);
        if ng_eq
            Gp_eq=g_eq(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);
            Gp_eq=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Gp_eq.^2));
        else
            Gp_eq=[];
        end
        Res=dX_quad-Fp;
        Res=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Res.^2));
        Res=scale_variables( [Res;Gp_eq]', data.dataNLP.data.discErrorConstScaling, 0 )';
        
        if mode==1
%             if ng_neq
%                 Gc=g(X,U,P,T,data.dataNLP.data);
%                 Gc_neq=Gc(:,ng_eq+1:end);
%             else
%                 Gc_neq=[];
%             end
%             g_vect=reshape(Gc_neq',M*ng,1);
            g_vect=reshape(g(X,U,P,T,data.dataNLP.data)',M*ng,1);
            if isfield(data.dataNLP.data,'Xscale')
                cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
            else
                cr=avrc(X,U,P,T,data.dataNLP)';
            end
            const=[ g_vect(data.dataNLP.gAllidx);
                    cr;
                    b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{});
                    sum(Res,2)];
        else
            const=[];
        end
end
      
end

