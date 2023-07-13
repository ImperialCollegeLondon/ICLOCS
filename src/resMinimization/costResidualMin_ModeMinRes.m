function [ varargout  ] = costResidualMin_ModeMinRes( X,U,P,T,data)
%costResidualMin_ModeMinRes - cost computation for integrated residual minimization (alternating method: residual minimization)
%
% Syntax:   [ ResNorm, Res_int ] = costResidualMin_ModeMinRes( X,U,P,T,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

% scale variable back
if isfield(data.dataNLP.data,'Xscale')
    X=scale_variables_back( X, data.dataNLP.data.Xscale_back, data.dataNLP.data.Xshift );
    U=scale_variables_back( U, data.dataNLP.data.Uscale_back, data.dataNLP.data.Ushift );
    if isfield(data.dataNLP.data,'Pscale')
        P=scale_variables_back( P, data.dataNLP.data.Pscale_back, data.dataNLP.data.Pshift );
    end
end

% Get function definitions
dataNLP=data.dataNLP;
f=dataNLP.data.InternalDynamics;

t0=T(1);tf=T(end);
if strcmp(dataNLP.options.discretization,'globalLGR') || strcmp(dataNLP.options.discretization,'hpLGR')
        % p/hp Transcription Method
        X_quad=(data.InterpH*X)./data.sumInterpH;
        U_quad=(data.InterpH*U)./data.sumInterpH;
        X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
        U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
        dX_quad=data.D_mat*X_quad./(data.DT_seg/2)/(tf-t0);
        
        P_quad=repelem(P(1,:),data.M_quad,1);
        
        ng_eq=dataNLP.sizes{18};
else
        % h Transcription Method
        n=dataNLP.sizes{3};

        F_k=f(X(1:2:end,:),U(1:2:end,:),P(1:2:end,:),data.tau(1:2:end)*(tf-t0),data.dataNLP.data);
        F_kph=data.DxHS_hf*X/(tf-t0)-F_k(1:end-1,:)/2;
        F_kp1=data.DxHS_p1*X/(tf-t0)+F_k(1:end-1,:);
        F=[F_k(1:end-1,:) F_kph F_kp1]';
        F=reshape(F(:),n,3*data.nps)';
        
        X_quad=repelem(X(1:2:end-1,:),data.npd_quad+1,1)+(tf-t0)*data.AxHS*F;
        U_quad=data.AuHS*U;
        dX_quad=data.AfHS*F;
        
        P_quad=repelem(P(1,:),data.M_quad,1);
        
        ng_eq=dataNLP.sizes{15};
end


Fp=f(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);        
Res=(dX_quad-Fp).^2;
Res_int=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Res));

if ng_eq
    Gp_eq=g_eq(X_quad,U_quad,P_quad,data.tau_quad*(tf-t0),data.dataNLP.data);
    Gp_eq=transpose((tf-t0)*data.DT_seg_node_mat./2*data.sum_nps_quad*(Gp_eq.^2));
else
    Gp_eq=[];
end

ResNorm=sum([Res_int;Gp_eq].*data.ResNormScale',1);

if nargout==1
    varargout{1} = ResNorm;
else
    Res_int=scale_variables( [Res_int;Gp_eq]', data.dataNLP.data.discErrorConstScaling, 0 )';
    varargout{1} = ResNorm;
    varargout{2} = Res_int;
end

end

