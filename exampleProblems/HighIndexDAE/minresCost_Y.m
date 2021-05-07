% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function [ResNorm_intsum,Res_intsum] = minresCost_Y(X,U,p,T,data)
global ADiGator_minresCost_Y
if isempty(ADiGator_minresCost_Y); ADiGator_LoadData(); end
Gator1Data = ADiGator_minresCost_Y.minresCost_Y.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %costResidualMin_ModeMinRes_Adigator - cost computation for integrated residual minimization (alternating method: residual minimization) with Adigator
%User Line: %
%User Line: % Syntax:   [ ResNorm_intsum, Res_intsum ] = costResidualMin_ModeMinRes_Adigator( X,U,p,T,data)
%User Line: %
%User Line: % Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the MIT License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.5
%User Line: % 1 Aug 2019
%User Line: % iclocs@imperial.ac.uk
dataNLP=data.dataNLP;
%User Line: dataNLP=data.dataNLP;
f=dataNLP.data.InternalDynamics;
%User Line: f=dataNLP.data.InternalDynamics;
dyn_data=data.dataNLP.data;
%User Line: dyn_data=data.dataNLP.data;
cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
%User Line: cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
    %User Line: % p/hp Transcription Method
    n.f = dataNLP.sizes{3};
    %User Line: n=dataNLP.sizes{3};
    ng_eq.f = dataNLP.sizes{18};
    %User Line: ng_eq=dataNLP.sizes{18};
    t_0.dY = T.dY(1);
    t_0.f = T.f(1);
    %User Line: t_0=T(1);
    cada1f1 = length(T.f);
    t_f.dY = T.dY(2);
    t_f.f = T.f(cada1f1);
    %User Line: t_f=T(end);
    cada1td1 = zeros(2,1);
    cada1td1(2) = t_f.dY;
    cada1td1(1) = cada1td1(1) + -t_0.dY;
    delta_t.dY = cada1td1;
    delta_t.f = t_f.f - t_0.f;
    %User Line: delta_t=t_f-t_0;
    cada1f1 = size(U.f,1);
    cada1f2dY = U.dY(Gator1Data.Index1);
    cada1f2 = U.f(cada1f1,:);
    cada1td1 = zeros(82,1);
    cada1td1(Gator1Data.Index2) = U.dY;
    cada1td1(Gator1Data.Index3) = cada1f2dY;
    U.dY = cada1td1;
    U.f = [U.f;cada1f2];
    %User Line: U=[U;U(end,:)];
    cada1td1 = sparse(Gator1Data.Index4,Gator1Data.Index5,X.dY,41,164);
    cada1td1 = data.InterpH*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f1dY = full(cada1td1(Gator1Data.Index6));
    cada1f1 = data.InterpH*X.f;
    cada1td1 = sparse(Gator1Data.Index7,Gator1Data.Index8,cada1f1dY,176,164);
    cada1td1 = data.sumInterpHMat*cada1td1;
    cada1td1 = cada1td1(:);
    X_quad.dY = full(cada1td1(Gator1Data.Index9));
    X_quad.f = data.sumInterpHMat*cada1f1;
    %User Line: X_quad=data.sumInterpHMat*(data.InterpH*X);
    cada1td1 = sparse(Gator1Data.Index10,Gator1Data.Index11,U.dY,41,80);
    cada1td1 = data.InterpH*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f1dY = full(cada1td1(Gator1Data.Index12));
    cada1f1 = data.InterpH*U.f;
    cada1td1 = sparse(Gator1Data.Index13,Gator1Data.Index14,cada1f1dY,176,80);
    cada1td1 = data.sumInterpHMat*cada1td1;
    cada1td1 = cada1td1(:);
    U_quad.dY = full(cada1td1(Gator1Data.Index15));
    U_quad.f = data.sumInterpHMat*cada1f1;
    %User Line: U_quad=data.sumInterpHMat*(data.InterpH*U);
    cada1f1dY = X.dY(Gator1Data.Index16);
    cada1f1 = X.f(data.interp_fixj,:);
    cada1td1 = zeros(3904,1);
    cada1td1(Gator1Data.Index17) = cada1f1dY;
    cada1td1(Gator1Data.Index18) = X_quad.dY(Gator1Data.Index19);
    X_quad.dY = cada1td1;
    X_quad.f(data.interp_fixi,:) = cada1f1;
    %User Line: X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
    cada1f1dY = U.dY(Gator1Data.Index20);
    cada1f1 = U.f(data.interp_fixj,:);
    cada1td1 = zeros(1912,1);
    cada1td1(Gator1Data.Index21) = cada1f1dY;
    cada1td1(Gator1Data.Index22) = U_quad.dY(Gator1Data.Index23);
    U_quad.dY = cada1td1;
    U_quad.f(data.interp_fixi,:) = cada1f1;
    %User Line: U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
    cada1td1 = sparse(Gator1Data.Index24,Gator1Data.Index25,X_quad.dY,176,164);
    cada1td1 = data.D_mat*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f1dY = full(cada1td1(Gator1Data.Index26));
    cada1f1 = data.D_mat*X_quad.f;
    cada1td1 = sparse(Gator1Data.Index27,Gator1Data.Index28,cada1f1dY,176,164);
    cada1td1 = data.DT_seg_mat_d2*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f2dY = full(cada1td1(Gator1Data.Index29));
    cada1f2 = data.DT_seg_mat_d2*cada1f1;
    cada1tempdY = delta_t.dY(Gator1Data.Index30);
    cada1td1 = zeros(5632,1);
    cada1td1(Gator1Data.Index31) = cada1f2dY./delta_t.f;
    cada1tf1 = cada1f2(Gator1Data.Index32);
    cada1td1(Gator1Data.Index33) = cada1td1(Gator1Data.Index33) + -cada1tf1(:)./delta_t.f.^2.*cada1tempdY;
    dX_quad.dY = cada1td1;
    dX_quad.f = cada1f2/delta_t.f;
    %User Line: dX_quad=data.DT_seg_mat_d2*(data.D_mat*X_quad)/delta_t;
    P_quad.f = repmat(p,176,1);
    %User Line: P_quad=repmat(p,data.M_quad,1);
    cada1tempdY = delta_t.dY(Gator1Data.Index34);
    cada1tf1 = data.tau_quad(Gator1Data.Index36);
    T_quad.dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index35);
    T_quad.f = data.tau_quad*delta_t.f;
    %User Line: T_quad=data.tau_quad*delta_t;
    cada1f1dY = X_quad.dY(Gator1Data.Index37);
    cada1f1 = X_quad.f(:,1);
    cada1f2dY = X_quad.dY(Gator1Data.Index38);
    cada1f2 = X_quad.f(:,2);
    cada1f3dY = X_quad.dY(Gator1Data.Index39);
    cada1f3 = X_quad.f(:,3);
    cada1f4dY = X_quad.dY(Gator1Data.Index40);
    cada1f4 = X_quad.f(:,4);
    cada1f5dY = U_quad.dY(Gator1Data.Index41);
    cada1f5 = U_quad.f(:,1);
    cada1f6dY = U_quad.dY(Gator1Data.Index42);
    cada1f6 = U_quad.f(:,2);
    cada1f7dY = -cada1f5dY;
    cada1f7 = uminus(cada1f5);
    cada1tf1 = cada1f1(Gator1Data.Index43);
    cada1td1 = zeros(1932,1);
    cada1td1(Gator1Data.Index44) = cada1tf1(:).*cada1f7dY;
    cada1tf1 = cada1f7(Gator1Data.Index45);
    cada1td1(Gator1Data.Index46) = cada1td1(Gator1Data.Index46) + cada1tf1(:).*cada1f1dY;
    cada1f8dY = cada1td1;
    cada1f8 = cada1f7.*cada1f1;
    cada1f9dY = dyn_data.a.*cada1f2dY;
    cada1f9 = dyn_data.a*cada1f2;
    cada1td1 = zeros(2908,1);
    cada1td1(Gator1Data.Index47) = cada1f8dY;
    cada1td1(Gator1Data.Index48) = cada1td1(Gator1Data.Index48) + -cada1f9dY;
    cada1f10dY = cada1td1;
    cada1f10 = cada1f8 - cada1f9;
    cada1tf1 = cada1f3(Gator1Data.Index49);
    cada1td1 = zeros(1932,1);
    cada1td1(Gator1Data.Index50) = cada1tf1(:).*cada1f6dY;
    cada1tf1 = cada1f6(Gator1Data.Index51);
    cada1td1(Gator1Data.Index52) = cada1td1(Gator1Data.Index52) + cada1tf1(:).*cada1f3dY;
    cada1f11dY = cada1td1;
    cada1f11 = cada1f6.*cada1f3;
    cada1td1 = zeros(4840,1);
    cada1td1(Gator1Data.Index53) = cada1f10dY;
    cada1td1(Gator1Data.Index54) = cada1td1(Gator1Data.Index54) + cada1f11dY;
    cada1f12dY = cada1td1;
    cada1f12 = cada1f10 + cada1f11;
    cada1f13 = uminus(dyn_data.g);
    cada1f14dY = 2.*cada1f5dY;
    cada1f14 = 2*cada1f5;
    cada1tf1 = cada1f3(Gator1Data.Index55);
    cada1td1 = zeros(1932,1);
    cada1td1(Gator1Data.Index56) = cada1tf1(:).*cada1f14dY;
    cada1tf1 = cada1f14(Gator1Data.Index57);
    cada1td1(Gator1Data.Index58) = cada1td1(Gator1Data.Index58) + cada1tf1(:).*cada1f3dY;
    cada1f15dY = cada1td1;
    cada1f15 = cada1f14.*cada1f3;
    cada1f16dY = -cada1f15dY;
    cada1f16 = cada1f13 - cada1f15;
    cada1tf1 = cada1f1(Gator1Data.Index59);
    cada1td1 = zeros(1932,1);
    cada1td1(Gator1Data.Index60) = cada1tf1(:).*cada1f6dY;
    cada1tf1 = cada1f6(Gator1Data.Index61);
    cada1td1(Gator1Data.Index62) = cada1td1(Gator1Data.Index62) + cada1tf1(:).*cada1f1dY;
    cada1f17dY = cada1td1;
    cada1f17 = cada1f6.*cada1f1;
    cada1td1 = zeros(3864,1);
    cada1td1(Gator1Data.Index63) = cada1f16dY;
    cada1td1(Gator1Data.Index64) = cada1td1(Gator1Data.Index64) + -cada1f17dY;
    cada1f18dY = cada1td1;
    cada1f18 = cada1f16 - cada1f17;
    cada1f19dY = dyn_data.a.*cada1f4dY;
    cada1f19 = dyn_data.a*cada1f4;
    cada1td1 = zeros(4840,1);
    cada1td1(Gator1Data.Index65) = cada1f18dY;
    cada1td1(Gator1Data.Index66) = cada1td1(Gator1Data.Index66) + -cada1f19dY;
    cada1f20dY = cada1td1;
    cada1f20 = cada1f18 - cada1f19;
    cada1tf2 = cada1f1(Gator1Data.Index67);
    cada1f21dY = 2.*cada1tf2(:).^(2-1).*cada1f1dY;
    cada1f21 = cada1f1.^2;
    cada1tf2 = cada1f3(Gator1Data.Index68);
    cada1f22dY = 2.*cada1tf2(:).^(2-1).*cada1f3dY;
    cada1f22 = cada1f3.^2;
    cada1td1 = zeros(1952,1);
    cada1td1(Gator1Data.Index69) = cada1f21dY;
    cada1td1(Gator1Data.Index70) = cada1td1(Gator1Data.Index70) + cada1f22dY;
    cada1f23dY = cada1td1;
    cada1f23 = cada1f21 + cada1f22;
    Gp.dY = cada1f23dY;
    Gp.f = cada1f23 - dyn_data.L;
    cada1td1 = zeros(11632,1);
    cada1td1(Gator1Data.Index71) = cada1f2dY;
    cada1td1(Gator1Data.Index72) = cada1f12dY;
    cada1td1(Gator1Data.Index73) = cada1f4dY;
    cada1td1(Gator1Data.Index74) = cada1f20dY;
    Fp.dY = cada1td1;
    Fp.f = [cada1f2 cada1f12 cada1f4 cada1f20];
    %User Line: [Fp,Gp]=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
    cada1f1 = 1:ng_eq.f;
    Gp.dY = Gp.dY(Gator1Data.Index75);
    Gp.f = Gp.f(:,cada1f1);
    %User Line: Gp=Gp(:,1:ng_eq);
    cada1td1 = zeros(15312,1);
    cada1td1(Gator1Data.Index76) = dX_quad.dY;
    cada1td1(Gator1Data.Index77) = cada1td1(Gator1Data.Index77) + -Fp.dY;
    cada1f1dY = cada1td1;
    cada1f1 = dX_quad.f - Fp.f;
    cada1tf2 = cada1f1(Gator1Data.Index78);
    cada1f2dY = 2.*cada1tf2(:).^(2-1).*cada1f1dY;
    cada1f2 = cada1f1.^2;
    cada1tf2 = Gp.f(Gator1Data.Index79);
    cada1f3dY = 2.*cada1tf2(:).^(2-1).*Gp.dY;
    cada1f3 = Gp.f.^2;
    cada1td1 = zeros(17264,1);
    cada1td1(Gator1Data.Index80) = cada1f2dY;
    cada1td1(Gator1Data.Index81) = cada1f3dY;
    Res.dY = cada1td1;
    Res.f = [cada1f2 cada1f3];
    %User Line: Res=[(dX_quad-Fp).^2 Gp.^2];
    cada1tempdY = delta_t.dY(Gator1Data.Index82);
    cada1tf1 = data.DT_seg_node_mat(Gator1Data.Index84);
    cada1f1dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index83);
    cada1f1 = delta_t.f*data.DT_seg_node_mat;
    cada1f2dY = cada1f1dY./2;
    cada1f2 = cada1f1/2;
    cada1td1 = zeros(8,16);
    cada1td1(Gator1Data.Index85) = cada1f2dY;
    cada1td1 = data.sum_nps_quad.'*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f3dY = cada1td1(Gator1Data.Index86);
    cada1f3 = cada1f2*data.sum_nps_quad;
    cada1td2 = sparse(Gator1Data.Index87,Gator1Data.Index88,cada1f3dY,176,16);
    cada1td2 = Res.f.'*cada1td2;
    cada1td1 = zeros(844,1);
    cada1td1(Gator1Data.Index90) = cada1td2(Gator1Data.Index89);
    cada1td2 = sparse(Gator1Data.Index91,Gator1Data.Index92,Res.dY,176,660);
    cada1td2 = cada1f3*cada1td2;
    cada1td2 = cada1td2(:);
    cada1td1(Gator1Data.Index94) = cada1td1(Gator1Data.Index94) + cada1td2(Gator1Data.Index93);
    Res_int.dY = cada1td1;
    Res_int.f = cada1f3*Res.f;
    %User Line: Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
    %User Line: % h Transcription Method
    %User Line: n=dataNLP.sizes{3};
    %User Line: nt=dataNLP.sizes{1};
    %User Line: M=dataNLP.sizes{7};
    %User Line: ng_eq=dataNLP.sizes{15};
    %User Line: cadaconditional1 = nt;
        %User Line: t_0=T(1);
        %User Line: t_f=T(end);
        %User Line: delta_t=t_f-t_0;
        %User Line: P=repmat(p,M,1);
        %User Line: X_col=X(1:2:end,:);
        %User Line: U_col=U(1:2:end,:);
        %User Line: T_col=data.tau(1:2:end)*delta_t;
        %User Line: F_k=f(X_col,U_col,P,T_col,dyn_data);
        %User Line: F_kph=data.DxHS_hf*X/delta_t-F_k(1:end-1,:)/2;
        %User Line: F_kp1=data.DxHS_p1*X/delta_t+F_k(1:end-1,:);
        %User Line: F=[F_k(1:end-1,:) F_kph F_kp1]';
        %User Line: F=reshape(F(:),n,3*data.nps)';
        %User Line: X_quad=data.repXend_mat*X(1:2:end-1,:)+delta_t*data.AxHS*F;
        %User Line: U_quad=data.AuHS*U;
        %User Line: dX_quad=data.AfHS*F;
        %User Line: P_quad=repmat(p,data.M_quad,1);
        %User Line: T_quad=data.tau_quad*delta_t;
        %User Line: [Fp,Gp]=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
        %User Line: Gp=Gp(:,1:ng_eq);
        %User Line: Res=[(dX_quad-Fp).^2 Gp.^2];
        %User Line: Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
        %User Line: t0=dataNLP.t0;
        %User Line: tf=dataNLP.tf;
        %User Line: deltat=tf-t0;
        %User Line: P=repmat(p,M,1);
        %User Line: X_col=X(1:2:end,:);
        %User Line: U_col=U(1:2:end,:);
        %User Line: T_col=data.tau(1:2:end)*deltat;
        %User Line: F_k=f(X_col,U_col,P,T_col,dyn_data);
        %User Line: F_kph=data.DxHS_hf*X/deltat-F_k(1:end-1,:)/2;
        %User Line: F_kp1=data.DxHS_p1*X/deltat+F_k(1:end-1,:);
        %User Line: F=[F_k(1:end-1,:) F_kph F_kp1]';
        %User Line: F=reshape(F(:),n,3*data.nps)';
        %User Line: X_quad=data.repXend_mat*X(1:2:end-1,:)+deltat*data.AxHS*F;
        %User Line: U_quad=data.AuHS*U;
        %User Line: dX_quad=data.AfHS*F;
        %User Line: P_quad=repmat(p,data.M_quad,1);
        %User Line: T_quad=data.tau_quad*deltat;
        %User Line: [Fp,Gp]=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
        %User Line: Gp=Gp(:,1:ng_eq);
        %User Line: Res=[(dX_quad-Fp).^2 Gp.^2];
        %User Line: Res_int=deltat*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
%User Line: % Compuation of integrated residual for each dynamics equation
cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1 = sparse(Gator1Data.Index95,Gator1Data.Index96,cada1f1dY,40,246);
cada1td1 = data.ResConstScaleMat*cada1td1;
cada1td1 = cada1td1(:);
Res_int_Const.dY = full(cada1td1(Gator1Data.Index97));
Res_int_Const.f = data.ResConstScaleMat*cada1f1;
%User Line: Res_int_Const=data.ResConstScaleMat*Res_int(:);
cada1f1 = n.f + ng_eq.f;
Res_int_Const.dY = Res_int_Const.dY;
Res_int_Const.f = reshape(Res_int_Const.f,data.nps,cada1f1);
%User Line: Res_int_Const=reshape(Res_int_Const,data.nps,n+ng_eq);
cada1td1 = sum(sparse(Gator1Data.Index98,Gator1Data.Index99,Res_int_Const.dY,8,662),1);
cada1f1dY = full(cada1td1(:));
cada1f1 = sum(Res_int_Const.f);
Res_intsum.dY = cada1f1dY;
Res_intsum.f = cada1f1.';
%User Line: Res_intsum=sum(Res_int_Const)';
%User Line: % Compuation of integrated residual norm
cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1 = sparse(Gator1Data.Index100,Gator1Data.Index101,cada1f1dY,40,246);
cada1td1 = data.ResNormScaleMat*cada1td1;
cada1td1 = cada1td1(:);
Res_int_Norm.dY = full(cada1td1(Gator1Data.Index102));
Res_int_Norm.f = data.ResNormScaleMat*cada1f1;
%User Line: Res_int_Norm=data.ResNormScaleMat*Res_int(:);
cada1f1 = n.f + ng_eq.f;
Res_int_Norm.dY = Res_int_Norm.dY;
Res_int_Norm.f = reshape(Res_int_Norm.f,data.nps,cada1f1);
%User Line: Res_int_Norm=reshape(Res_int_Norm,data.nps,n+ng_eq);
cada1td1 = sum(sparse(Gator1Data.Index103,Gator1Data.Index104,Res_int_Norm.dY,8,662),1);
cada1f1dY = full(cada1td1(:));
cada1f1 = sum(Res_int_Norm.f);
cada1td1 = sum(sparse(Gator1Data.Index105,Gator1Data.Index106,cada1f1dY,5,246),1);
ResNorm_intsum.dY = full(cada1td1(:));
ResNorm_intsum.f = sum(cada1f1);
%User Line: ResNorm_intsum=sum(sum(Res_int_Norm));
ResNorm_intsum.dY_size = 246;
ResNorm_intsum.dY_location = Gator1Data.Index107;
Res_intsum.dY_size = [5,246];
Res_intsum.dY_location = Gator1Data.Index108;
end


function ADiGator_LoadData()
global ADiGator_minresCost_Y
ADiGator_minresCost_Y = load('minresCost_Y.mat');
return
end