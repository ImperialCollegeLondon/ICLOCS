function [ const ] = constResidualMin_ModeMinRes( X, U, P, T, data)
%constResidualMin_ModeMinRes - constraint function computation for
%integrated residual minimization (alternating method: residual
%minimization)
%
% Syntax:   [ const ] = constResidualMin_ModeMinRes( X, U, P, T, data)
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

if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
    X_endpoint=X([1,end],:);
    x_endpoint=data.dataNLP.data.x_endpoint([1,end],:);
end

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

x0=X(1,:);
u0=U(1,:)';
xf=X(end,:)';
uf=U(end,:)';

[L,E,f,g,avrc,b]=deal(data.dataNLP.functions_unscaled{:});

switch dataNLP.options.discretization
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        [nt,np,n,m,ng,~,M,~,~,~,~,~,~,~,~,~,~]=deal(dataNLP.sizes{1:17});
        X_mesh=X(1:end-1,:);
        U_mesh=U(1:end-1,:);
        T_mesh=T(1:end-1,:);
        X_mesh_Np1=X;
        
        g_vect=reshape(g(X_mesh,U_mesh,P,T_mesh,data.dataNLP.data)',M*ng,1);
        Cost=data.dataNLP.map.w'*((tf-t0)*data.DT_seg_mesh./2.*L(X_mesh,Xr,U_mesh,Ur,P,T_mesh,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
        if isfield(dataNLP.options,'resminRep') && isfield(dataNLP.options.resminRep,'costTol')
            costTol=dataNLP.options.resminRep.costTol;
        else
            costTol=0;
        end
        if isfield(data,'cost_org')
            if data.cost_org<0
                Costconst= Cost*(1+costTol);
            else
                Costconst= Cost*(1-costTol);
            end
        else
            Costconst= Cost;
        end
        
        if isfield(data.dataNLP.data,'Xscale')
            cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
        else
            cr=avrc(X_mesh_Np1,U_mesh,P,T,data.dataNLP)';
        end


        if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
            const=[  g_vect(data.dataNLP.gAllidx);
                     cr;
                     X_endpoint(:)-x_endpoint(:);
                     Costconst];
        else
            const=[  g_vect(data.dataNLP.gAllidx);
                     cr;
                     b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{});
                     Costconst];
        end
    otherwise
        [nt,np,n,m,ng,~,M,N,ns,~,~,~,~]=deal(dataNLP.sizes{1:13});
        g_vect=reshape(g(X,U,P,T,data.dataNLP.data)',M*ng,1);
        Cost=data.dataNLP.map.w'*((tf-t0)*L(X,[],U,[],P,T,data.dataNLP.data))+E(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data);
        if isfield(dataNLP.options,'resminRep') && isfield(dataNLP.options.resminRep,'costTol')
            costTol=dataNLP.options.resminRep.costTol;
        else
            costTol=0;
        end
        if isfield(data,'cost_org')
            if data.cost_org<0
                Costconst= Cost*(1+costTol);
            else
                Costconst= Cost*(1-costTol);
            end
        else
            Costconst= Cost;
        end
        
        if isfield(data.dataNLP.data,'Xscale')
            cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
        else
            cr=avrc(X,U,P,T,data.dataNLP)';
        end
        
        if isfield(dataNLP.options,'resminRep') && dataNLP.options.resminRep.collmatch
            const=[ g_vect(data.dataNLP.gAllidx);
                    cr;
                    X_endpoint(:)-x_endpoint(:);
                    Costconst];
        else
            const=[ g_vect(data.dataNLP.gAllidx);
                    cr;
                    b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{});
                    Costconst];
        end
end
      
end

