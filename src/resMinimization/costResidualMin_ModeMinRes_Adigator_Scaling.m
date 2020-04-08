function [ ResNorm_intsum, Res_intsum ] = costResidualMin_ModeMinRes_Adigator_Scaling( X,U,p,T,data)
%costResidualMin_ModeMinRes_Adigator - cost computation for integrated
%residual minimization (alternating method: residual minimization) with
%Adigator and automatic scaling
%
% Syntax:   [ ResNorm_intsum, Res_intsum ] = costResidualMin_ModeMinRes_Adigator_Scaling( X,U,p,T,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


dataNLP=data.dataNLP;
f=dataNLP.data.InternalDynamics;
dyn_data=data.dataNLP.data;

% unscaling of variables
dimX1=size(X,1);dimX2=size(X,2);
dimU1=size(U,1);dimU2=size(U,2);
XshiftMat=dataNLP.scaling.XshiftMat;
UshiftMat=dataNLP.scaling.UshiftMat;
X=dataNLP.scaling.XunscaleMat*(X(:)-XshiftMat(:));
U=dataNLP.scaling.UunscaleMat*(U(:)-UshiftMat(:));
X=reshape(X,dimX1,dimX2);
U=reshape(U,dimU1,dimU2);
if isfield(dataNLP.data,'Pscale')
    p=(p-dataNLP.data.Pshift)./dataNLP.data.Pscale;
end

        
if strcmp(dataNLP.options.discretization,'globalLGR') || strcmp(dataNLP.options.discretization,'hpLGR')
	% p/hp Transcription Method
        n=dataNLP.sizes{3};
        
        t_0=T(1);
        t_f=T(end);
        delta_t=t_f-t_0;
        U=[U;U(end,:)];
        X_quad=data.sumInterpHMat*(data.InterpH*X);
        U_quad=data.sumInterpHMat*(data.InterpH*U);
        X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
        U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
        dX_quad=data.DT_seg_mat_d2*(data.D_mat*X_quad)/delta_t;
        

        P_quad=repmat(p,data.M_quad,1);
        T_quad=data.tau_quad*delta_t;
        Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
        Res=(dX_quad-Fp).^2;
        Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
else
	% h Transcription Method
        n=dataNLP.sizes{3};
        nt=dataNLP.sizes{1};
        M=dataNLP.sizes{7};
        if nt % free time
            t_0=T(1);
            t_f=T(end);
            delta_t=t_f-t_0;
            
            P=repmat(p,M,1);
            X_col=X(1:2:end,:);
            U_col=U(1:2:end,:);
            T_col=data.tau(1:2:end)*delta_t;


            F_k=f(X_col,U_col,P,T_col,dyn_data);
            F_kph=data.DxHS_hf*X/delta_t-F_k(1:end-1,:)/2;
            F_kp1=data.DxHS_p1*X/delta_t+F_k(1:end-1,:);
            F=[F_k(1:end-1,:) F_kph F_kp1]';
            F=reshape(F(:),n,3*data.nps)';

            X_quad=data.repXend_mat*X(1:2:end-1,:)+delta_t*data.AxHS*F;
            U_quad=data.AuHS*U;
            dX_quad=data.AfHS*F;

            P_quad=repmat(p,data.M_quad,1);
            T_quad=data.tau_quad*delta_t;
            
            Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
            Res=(dX_quad-Fp).^2;
            Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
        else % fixed time
            t0=dataNLP.t0;
            tf=dataNLP.tf;
            deltat=tf-t0;
            
            P=repmat(p,M,1);
            X_col=X(1:2:end,:);
            U_col=U(1:2:end,:);
            T_col=data.tau(1:2:end)*deltat;


            F_k=f(X_col,U_col,P,T_col,dyn_data);
            F_kph=data.DxHS_hf*X/deltat-F_k(1:end-1,:)/2;
            F_kp1=data.DxHS_p1*X/deltat+F_k(1:end-1,:);
            F=[F_k(1:end-1,:) F_kph F_kp1]';
            F=reshape(F(:),n,3*data.nps)';

            X_quad=data.repXend_mat*X(1:2:end-1,:)+deltat*data.AxHS*F;
            U_quad=data.AuHS*U;
            dX_quad=data.AfHS*F;

            P_quad=repmat(p,data.M_quad,1);
            T_quad=data.tau_quad*deltat;
            
            Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
            Res=(dX_quad-Fp).^2;
            Res_int=deltat*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;

        end

end

% Compuation of integrated residual for each dynamics equation
Res_int_Const=data.ResConstScaleMat*Res_int(:);
Res_int_Const=reshape(Res_int_Const,data.nps,n);
Res_intsum=sum(Res_int_Const)';

% Compuation of integrated residual norm
Res_int_Norm=data.ResNormScaleMat*Res_int(:);
Res_int_Norm=reshape(Res_int_Norm,data.nps,n);
ResNorm_intsum=sum(sum(Res_int_Norm));


end

