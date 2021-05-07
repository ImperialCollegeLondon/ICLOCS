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

function [ResNorm_intsum,Res_intsum] = minresCost_YY(X,U,p,T,data)
global ADiGator_minresCost_YY
if isempty(ADiGator_minresCost_YY); ADiGator_LoadData(); end
Gator1Data = ADiGator_minresCost_YY.minresCost_YY.Gator1Data;
Gator2Data = ADiGator_minresCost_YY.minresCost_YY.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: %costResidualMin_ModeMinRes_Adigator - cost computation for integrated
%User Line: %residual minimization (alternating method: residual minimization) with
%User Line: %Adigator and automatic scaling
%User Line: %
%User Line: % Syntax:   [ ResNorm_intsum, Res_intsum ] = costResidualMin_ModeMinRes_Adigator_Scaling( X,U,p,T,data)
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
%User Line: % unscaling of variables
dimX1.f = size(X.f,1);
%User Line: dimX1=size(X,1);
dimX2.f = size(X.f,2);
%User Line: dimX2=size(X,2);
dimU1.f = size(U.f,1);
%User Line: dimU1=size(U,1);
dimU2.f = size(U.f,2);
%User Line: dimU2=size(U,2);
XshiftMat = dataNLP.scaling.XshiftMat;
%User Line: XshiftMat=dataNLP.scaling.XshiftMat;
UshiftMat = dataNLP.scaling.UshiftMat;
%User Line: UshiftMat=dataNLP.scaling.UshiftMat;
%User Line: % X=dataNLP.scaling.XunscaleMat*(X(:)-XshiftMat(:));
%User Line: % U=dataNLP.scaling.UunscaleMat*(U(:)-UshiftMat(:));
cada1f1dY = X.dY;
cada1f1 = X.f(:);
cada1td1 = sparse(Gator1Data.Index1,Gator1Data.Index2,cada1f1dY,597,597);
cada1td1 = dataNLP.scaling.XunscaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index3);
cada1f2dY = full(cada2f1);
cada1f2 = dataNLP.scaling.XunscaleMat*cada1f1;
cada1f3 = XshiftMat(:);
X.dY = cada1f2dY;
X.f = cada1f2 - cada1f3;
%User Line: X=dataNLP.scaling.XunscaleMat*X(:)-XshiftMat(:);
cada1f1dY = U.dY;
cada1f1 = U.f(:);
cada1td1 = sparse(Gator1Data.Index4,Gator1Data.Index5,cada1f1dY,199,199);
cada1td1 = dataNLP.scaling.UunscaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index6);
cada1f2dY = full(cada2f1);
cada1f2 = dataNLP.scaling.UunscaleMat*cada1f1;
cada1f3 = UshiftMat(:);
U.dY = cada1f2dY;
U.f = cada1f2 - cada1f3;
%User Line: U=dataNLP.scaling.UunscaleMat*U(:)-UshiftMat(:);
X.dY = X.dY;
X.f = reshape(X.f,dimX1.f,dimX2.f);
%User Line: X=reshape(X,dimX1,dimX2);
U.dY = U.dY;
U.f = reshape(U.f,dimU1.f,dimU2.f);
%User Line: U=reshape(U,dimU1,dimU2);
cadaconditional1 = isfield(dataNLP.data,'Pscale');
%User Line: cadaconditional1 = isfield(dataNLP.data,'Pscale');
%User Line: %     p=(p-dataNLP.data.Pshift)./dataNLP.data.Pscale;
%User Line: p=p./dataNLP.data.Pscale-dataNLP.data.Pshift;
cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
%User Line: cadaconditional1 = strcmp(dataNLP.options.discretization,'globalLGR') | strcmp(dataNLP.options.discretization,'hpLGR');
%User Line: % p/hp Transcription Method
%User Line: n=dataNLP.sizes{3};
%User Line: t_0=T(1);
%User Line: t_f=T(end);
%User Line: delta_t=t_f-t_0;
%User Line: U=[U;U(end,:)];
%User Line: X_quad=data.sumInterpHMat*(data.InterpH*X);
%User Line: U_quad=data.sumInterpHMat*(data.InterpH*U);
%User Line: X_quad(data.interp_fixi,:)=X(data.interp_fixj,:);
%User Line: U_quad(data.interp_fixi,:)=U(data.interp_fixj,:);
%User Line: dX_quad=data.DT_seg_mat_d2*(data.D_mat*X_quad)/delta_t;
%User Line: P_quad=repmat(p,data.M_quad,1);
%User Line: T_quad=data.tau_quad*delta_t;
%User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
%User Line: Res=(dX_quad-Fp).^2;
%User Line: Res_int=delta_t*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
%User Line: % h Transcription Method
n.f = dataNLP.sizes{3};
%User Line: n=dataNLP.sizes{3};
nt.f = dataNLP.sizes{1};
%User Line: nt=dataNLP.sizes{1};
M.f = dataNLP.sizes{7};
%User Line: M=dataNLP.sizes{7};
cadaconditional1 = nt.f;
%User Line: cadaconditional1 = nt;
t_0.dY = T.dY(1);
t_0.f = T.f(1);
%User Line: t_0=T(1);
cada1f1 = length(T.f);
t_f.dY = T.dY(2);
t_f.f = T.f(cada1f1);
%User Line: t_f=T(end);
cada1td1 =  zeros(2,1);
cada1td1(2) = t_f.dY;
cada2f1 = cada1td1(1);
cada2f2 = uminus(t_0.dY);
cada2f3 = cada2f1 + cada2f2;
cada1td1(1) = cada2f3;
delta_t.dY = cada1td1;
delta_t.f = t_f.f - t_0.f;
%User Line: delta_t=t_f-t_0;
P.f = repmat(p,199,1);
%User Line: P=repmat(p,M,1);
cada1f1 = size(X.f,1);
cada1f2 = 1:2:cada1f1;
X_col.dY = X.dY(Gator1Data.Index7);
X_col.f = X.f(cada1f2,:);
%User Line: X_col=X(1:2:end,:);
cada1f1 = size(U.f,1);
cada1f2 = 1:2:cada1f1;
U_col.dY = U.dY(Gator1Data.Index8);
U_col.f = U.f(cada1f2,:);
%User Line: U_col=U(1:2:end,:);
cada1f1 = length(data.tau);
cada1f2 = 1:2:cada1f1;
cada1f3 = data.tau(cada1f2);
cada1tempdY = delta_t.dY(Gator1Data.Index9);
cada1tf1 = cada1f3(Gator1Data.Index11);
cada2f1 = cada1tf1(:);
cada2f2 = cada1tempdY(Gator1Data.Index10);
T_col.dY = cada2f1.*cada2f2;
T_col.f = cada1f3*delta_t.f;
%User Line: T_col=data.tau(1:2:end)*delta_t;
cada1f1dY = X_col.dY(Gator1Data.Index12);
cada1f1 = X_col.f(:,1);
cada1f2dY = X_col.dY(Gator1Data.Index13);
cada1f2 = X_col.f(:,2);
cada1f3dY = X_col.dY(Gator1Data.Index14);
cada1f3 = X_col.f(:,3);
cada1f4dY = U_col.dY(Gator1Data.Index15);
cada1f4 = U_col.f(:,1);
cada2f1dY = cada1f3dY;
cada2f1 = cada1f3(:);
cada2f2dY = 2.*cada2f1(:).^(2-1).*cada2f1dY;
cada2f2 = cada2f1.^2;
cada2f3dY = -(-1)./cada2f2(:).^2.*cada2f2dY;
cada2f3 = (-1)./cada2f2;
cada1f5dYdY = cada1f3dY(:).*cada2f3dY;
cada1f5dY = cada2f3.*cada1f3dY;
cada1f5 = 1./cada1f3;
cada2f1dY = cada1f2dY;
cada2f1 = cada1f2(:);
cada2f2dY = 1.*cada2f1(:).^(1-1).*cada2f1dY;
cada2f2 = cada2f1.^1;
cada2f3dY = 2.*cada2f2dY;
cada2f3 = 2*cada2f2;
cada1f6dYdY = cada1f2dY(:).*cada2f3dY;
cada1f6dY = cada2f3.*cada1f2dY;
cada1f6 = cada1f2.^2;
cada1f7dYdY = dyn_data.D0.*cada1f6dYdY;
cada1f7dY = dyn_data.D0*cada1f6dY;
cada1f7 = dyn_data.D0*cada1f6;
cada1f8dY = uminus(cada1f1dY);
cada1f8 = uminus(cada1f1);
cada1f9dY = cada1f8dY/dyn_data.H;
cada1f9 = cada1f8/dyn_data.H;
cada2f1dY = cada1f9dY;
cada2f1 = cada1f9(:);
cada2f2dY = exp(cada2f1(:)).*cada2f1dY;
cada2f2 = exp(cada2f1);
cada1f10dYdY = cada1f9dY(:).*cada2f2dY;
cada1f10dY = cada2f2.*cada1f9dY;
cada1f10 = exp(cada1f9);
cada1td1 =  zeros(200,1);
cada2f1dY = cada1f10dY;
cada2f1 = cada1f10(:);
cada2td1 = zeros(200,1);
cada2td1(Gator2Data.Index1) = cada1f7dY(:).*cada2f1dY;
cada2td1(Gator2Data.Index2) = cada2td1(Gator2Data.Index2) + cada2f1(:).*cada1f7dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f7dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index16) = cada2f2;
cada2f1 = cada1td1(Gator1Data.Index17);
cada2f2dY = cada1f7dY;
cada2f2 = cada1f7(:);
cada2td1 = zeros(200,1);
cada2td1(Gator2Data.Index3) = cada1f10dY(:).*cada2f2dY;
cada2td1(Gator2Data.Index4) = cada2td1(Gator2Data.Index4) + cada2f2(:).*cada1f10dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f10dY;
cada2f4dY = cada2f3dY;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(400,1);
cada2td1(Gator2Data.Index5) = cada2f4dY;
cada2td1(Gator2Data.Index6) = cada1td1dY(Gator2Data.Index7);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index17) = cada2f4;
cada1f11dYdY = cada1td1dY; cada1f11dY = cada1td1;
cada1f11 = cada1f7.*cada1f10;
cada1td1 =  zeros(300,1);
cada1td1(Gator1Data.Index18) = cada1f4dY;
cada2f1 = cada1td1(Gator1Data.Index19);
cada2f2dY = -cada1f11dYdY;
cada2f2 = uminus(cada1f11dY);
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada1td1dY = cada2f3dY;
cada1td1(Gator1Data.Index19) = cada2f3;
cada1f12dYdY = cada1td1dY; cada1f12dY = cada1td1;
cada1f12 = cada1f4 - cada1f11;
cada1td1 =  zeros(400,1);
cada2f1dY = cada1f12dY;
cada2f1 = cada1f12(:);
cada2tf1 = cada1f5dY(Gator2Data.Index8);
cada2td1 = zeros(400,1);
cada2td1(Gator2Data.Index9) = cada2tf1(:).*cada2f1dY;
cada2td1(Gator2Data.Index10) = cada2td1(Gator2Data.Index10) + cada2f1(:).*cada1f5dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f5dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index20) = cada2f2;
cada1tf1dY = cada1f5dY(Gator2Data.Index11);
cada1tf1 = cada1f5(Gator1Data.Index21);
cada2f1 = cada1td1(Gator1Data.Index22);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2td1 = zeros(700,1);
cada2td1(Gator2Data.Index12) = cada1f12dY(:).*cada2f2dY;
cada2tf1 = cada2f2(Gator2Data.Index13);
cada2td1(Gator2Data.Index14) = cada2td1(Gator2Data.Index14) + cada2tf1(:).*cada1f12dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f12dY;
cada2f4dY = cada2f3dY;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(1100,1);
cada2td1(Gator2Data.Index15) = cada2f4dY;
cada2td1(Gator2Data.Index16) = cada1td1dY(Gator2Data.Index17);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index22) = cada2f4;
cada1f13dYdY = cada1td1dY; cada1f13dY = cada1td1;
cada1f13 = cada1f5.*cada1f12;
cada1f14dYdY = cada1f13dYdY; cada1f14dY = cada1f13dY;
cada1f14 = cada1f13 - dyn_data.grav;
cada1f15dY = uminus(cada1f4dY);
cada1f15 = uminus(cada1f4);
cada1f16dY = cada1f15dY/dyn_data.c;
cada1f16 = cada1f15/dyn_data.c;
cada1td1 =  zeros(600,1);
cada1td1(Gator1Data.Index23) = cada1f2dY;
cada1td1dY = cada1f14dYdY;
cada1td1(Gator1Data.Index24) = cada1f14dY;
cada1td1(Gator1Data.Index25) = cada1f16dY;
F_k.dYdY = cada1td1dY; F_k.dY = cada1td1;
F_k.f = [cada1f2 cada1f14 cada1f16];
%User Line: F_k=f(X_col,U_col,P,T_col,dyn_data);
cada1td1 = sparse(Gator1Data.Index26,Gator1Data.Index27,X.dY,199,597);
cada1td1 = data.DxHS_hf*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index28);
cada1f1dY = full(cada2f1);
cada1f1 = data.DxHS_hf*X.f;
cada1tempdY = delta_t.dY(Gator1Data.Index29);
cada1td1 =  zeros(1485,1);
cada2tempdY = delta_t.dY(Gator2Data.Index18);
cada2tf1 = cada1f1dY(Gator2Data.Index19);
cada2f1dY = -cada2tf1(:)./delta_t.f.^2.*cada2tempdY;
cada2f1 = cada1f1dY/delta_t.f;
cada1td1dY = cada2f1dY;
cada1td1(Gator1Data.Index30) = cada2f1;
cada1tf1dY = cada1f1dY(Gator2Data.Index20);
cada1tf1 = cada1f1(Gator1Data.Index31);
cada2f1 = cada1td1(Gator1Data.Index32);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2f3dY = -cada2f2dY;
cada2f3 = uminus(cada2f2);
cada2f4dY = 2.*delta_t.f.^(2-1).*delta_t.dY;
cada2f4 = delta_t.f^2;
cada2tempdY = cada2f4dY(Gator2Data.Index21);
cada2td1 = zeros(2970,1);
cada2td1(Gator2Data.Index22) = cada2f3dY./cada2f4;
cada2tf1 = cada2f3(Gator2Data.Index23);
cada2td1(Gator2Data.Index24) = cada2td1(Gator2Data.Index24) + -cada2tf1(:)./cada2f4.^2.*cada2tempdY;
cada2f5dY = cada2td1;
cada2f5 = cada2f3/cada2f4;
cada2tf1 = cada1tempdY(Gator2Data.Index25);
cada2f6dY = cada2tf1(:).*cada2f5dY;
cada2f6 = cada2f5.*cada1tempdY;
cada2f7dY = cada2f6dY;
cada2f7 = cada2f1 + cada2f6;
cada2td1 = zeros(4752,1);
cada2td1(Gator2Data.Index26) = cada2f7dY;
cada2td1(Gator2Data.Index27) = cada1td1dY(Gator2Data.Index28);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index32) = cada2f7;
cada1f2dYdY = cada1td1dY; cada1f2dY = cada1td1;
cada1f2 = cada1f1/delta_t.f;
cada1f3 = size(F_k.f,1);
cada1f4 = cada1f3 - 1;
cada1f5 = 1:cada1f4;
cada1f6dYdY = F_k.dYdY(Gator2Data.Index29);
cada1f6dY = F_k.dY(Gator1Data.Index33);
cada1f6 = F_k.f(cada1f5,:);
cada1f7dYdY = cada1f6dYdY./2;
cada1f7dY = cada1f6dY/2;
cada1f7 = cada1f6/2;
cada1td1 =  zeros(1980,1);
cada1td1dY = cada1f2dYdY;
cada1td1(Gator1Data.Index34) = cada1f2dY;
cada2f1dY = cada1td1dY(Gator2Data.Index30);
cada2f1 = cada1td1(Gator1Data.Index35);
cada2f2dY = -cada1f7dYdY;
cada2f2 = uminus(cada1f7dY);
cada2td1 = zeros(1287,1);
cada2td1(Gator2Data.Index31) = cada2f1dY;
cada2td1(Gator2Data.Index32) = cada2td1(Gator2Data.Index32) + cada2f2dY;
cada2f3dY = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(5841,1);
cada2td1(Gator2Data.Index33) = cada2f3dY;
cada2td1(Gator2Data.Index34) = cada1td1dY(Gator2Data.Index35);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index35) = cada2f3;
F_kph.dYdY = cada1td1dY; F_kph.dY = cada1td1;
F_kph.f = cada1f2 - cada1f7;
%User Line: F_kph=data.DxHS_hf*X/delta_t-F_k(1:end-1,:)/2;
cada1td1 = sparse(Gator1Data.Index36,Gator1Data.Index37,X.dY,199,597);
cada1td1 = data.DxHS_p1*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index38);
cada1f1dY = full(cada2f1);
cada1f1 = data.DxHS_p1*X.f;
cada1tempdY = delta_t.dY(Gator1Data.Index39);
cada1td1 =  zeros(1485,1);
cada2tempdY = delta_t.dY(Gator2Data.Index36);
cada2tf1 = cada1f1dY(Gator2Data.Index37);
cada2f1dY = -cada2tf1(:)./delta_t.f.^2.*cada2tempdY;
cada2f1 = cada1f1dY/delta_t.f;
cada1td1dY = cada2f1dY;
cada1td1(Gator1Data.Index40) = cada2f1;
cada1tf1dY = cada1f1dY(Gator2Data.Index38);
cada1tf1 = cada1f1(Gator1Data.Index41);
cada2f1 = cada1td1(Gator1Data.Index42);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2f3dY = -cada2f2dY;
cada2f3 = uminus(cada2f2);
cada2f4dY = 2.*delta_t.f.^(2-1).*delta_t.dY;
cada2f4 = delta_t.f^2;
cada2tempdY = cada2f4dY(Gator2Data.Index39);
cada2td1 = zeros(2970,1);
cada2td1(Gator2Data.Index40) = cada2f3dY./cada2f4;
cada2tf1 = cada2f3(Gator2Data.Index41);
cada2td1(Gator2Data.Index42) = cada2td1(Gator2Data.Index42) + -cada2tf1(:)./cada2f4.^2.*cada2tempdY;
cada2f5dY = cada2td1;
cada2f5 = cada2f3/cada2f4;
cada2tf1 = cada1tempdY(Gator2Data.Index43);
cada2f6dY = cada2tf1(:).*cada2f5dY;
cada2f6 = cada2f5.*cada1tempdY;
cada2f7dY = cada2f6dY;
cada2f7 = cada2f1 + cada2f6;
cada2td1 = zeros(4752,1);
cada2td1(Gator2Data.Index44) = cada2f7dY;
cada2td1(Gator2Data.Index45) = cada1td1dY(Gator2Data.Index46);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index42) = cada2f7;
cada1f2dYdY = cada1td1dY; cada1f2dY = cada1td1;
cada1f2 = cada1f1/delta_t.f;
cada1f3 = size(F_k.f,1);
cada1f4 = cada1f3 - 1;
cada1f5 = 1:cada1f4;
cada1f6dYdY = F_k.dYdY(Gator2Data.Index47);
cada1f6dY = F_k.dY(Gator1Data.Index43);
cada1f6 = F_k.f(cada1f5,:);
cada1td1 =  zeros(1980,1);
cada1td1dY = cada1f2dYdY;
cada1td1(Gator1Data.Index44) = cada1f2dY;
cada2f1dY = cada1td1dY(Gator2Data.Index48);
cada2f1 = cada1td1(Gator1Data.Index45);
cada2td1 = zeros(1287,1);
cada2td1(Gator2Data.Index49) = cada2f1dY;
cada2td1(Gator2Data.Index50) = cada2td1(Gator2Data.Index50) + cada1f6dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1 + cada1f6dY;
cada2td1 = zeros(5841,1);
cada2td1(Gator2Data.Index51) = cada2f2dY;
cada2td1(Gator2Data.Index52) = cada1td1dY(Gator2Data.Index53);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index45) = cada2f2;
F_kp1.dYdY = cada1td1dY; F_kp1.dY = cada1td1;
F_kp1.f = cada1f2 + cada1f6;
%User Line: F_kp1=data.DxHS_p1*X/delta_t+F_k(1:end-1,:);
cada1f1 = size(F_k.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
cada1f4dYdY = F_k.dYdY(Gator2Data.Index54);
cada1f4dY = F_k.dY(Gator1Data.Index46);
cada1f4 = F_k.f(cada1f3,:);
cada1td1 =  zeros(4554,1);
cada1td1dY = cada1f4dYdY;
cada1td1(Gator1Data.Index47) = cada1f4dY;
cada2td1 = zeros(6930,1);
cada2td1(Gator2Data.Index55) = F_kph.dYdY;
cada2td1(Gator2Data.Index56) = cada1td1dY(Gator2Data.Index57);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index48) = F_kph.dY;
cada2td1 = zeros(12771,1);
cada2td1(Gator2Data.Index58) = F_kp1.dYdY;
cada2td1(Gator2Data.Index59) = cada1td1dY(Gator2Data.Index60);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index49) = F_kp1.dY;
cada1f5dYdY = cada1td1dY; cada1f5dY = cada1td1;
cada1f5 = [cada1f4 F_kph.f F_kp1.f];
F.dYdY = cada1f5dYdY(Gator2Data.Index61);
F.dY = cada1f5dY(Gator1Data.Index50);
F.f = cada1f5.';
%User Line: F=[F_k(1:end-1,:) F_kph F_kp1]';
cada1f1dYdY = F.dYdY; cada1f1dY = F.dY;
cada1f1 = F.f(:);
cada1f2 = 3*data.nps;
cada1f3dYdY = cada1f1dYdY; cada1f3dY = cada1f1dY;
cada1f3 = reshape(cada1f1,n.f,cada1f2);
F.dYdY = cada1f3dYdY(Gator2Data.Index62);
F.dY = cada1f3dY(Gator1Data.Index51);
F.f = cada1f3.';
%User Line: F=reshape(F(:),n,3*data.nps)';
cada1f1 = size(X.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:2:cada1f2;
cada1f4dY = X.dY(Gator1Data.Index52);
cada1f4 = X.f(cada1f3,:);
cada1td1 = sparse(Gator1Data.Index53,Gator1Data.Index54,cada1f4dY,99,297);
cada1td1 = data.repXend_mat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index55);
cada1f5dY = full(cada2f1);
cada1f5 = data.repXend_mat*cada1f4;
cada1tempdY = delta_t.dY(Gator1Data.Index56);
cada1tf1 = data.AxHS(Gator1Data.Index58);
cada2f1 = cada1tf1(:);
cada2f2 = cada1tempdY(Gator1Data.Index57);
cada1f6dY = cada2f1.*cada2f2;
cada1f6 = delta_t.f*data.AxHS;
cada1td2 = sparse(Gator1Data.Index59,Gator1Data.Index60,cada1f6dY,297,2644);
cada2f1dY = F.dY(Gator2Data.Index63);
cada2f1 = F.f.';
cada2td1 = sparse(Gator2Data.Index64,Gator2Data.Index65,cada2f1dY,297,1098);
cada2td1 = cada1td2.'*cada2td1;
cada2td1 = cada2td1(:);
cada1td2dY = full(cada2td1(Gator2Data.Index66));
cada1td2 = cada2f1*cada1td2;
cada1td1 =  zeros(26440,1);
cada2f1dY = cada1td2dY(Gator2Data.Index67);
cada2f1 = cada1td2(Gator1Data.Index61);
cada1td1dY = cada2f1dY;
cada1td1(Gator1Data.Index62) = cada2f1;
cada1td2dY = F.dYdY(Gator2Data.Index68);
cada1td2 = sparse(Gator1Data.Index63,Gator1Data.Index64,F.dY,297,1098);
cada2td2 = sparse(Gator2Data.Index69,Gator2Data.Index70,cada1f6dY,297,2644);
cada2td2 = cada1td2.'*cada2td2;
cada2td1 = zeros(91218,1);
cada2td1(Gator2Data.Index72) = cada2td2(Gator2Data.Index71);
cada2td2 = sparse(Gator2Data.Index73,Gator2Data.Index74,cada1td2dY,297,3489);
cada2td2 = cada1f6*cada2td2;
cada2td2 = cada2td2(:);
cada2td1(Gator2Data.Index76) = cada2td1(Gator2Data.Index76) + cada2td2(Gator2Data.Index75);
cada1td2dY = cada2td1;
cada1td2 = cada1f6*cada1td2;
cada1td2 = cada1td2(:);
cada2f1dY = cada1td2dY(Gator2Data.Index77);
cada2f1 = cada1td2(Gator1Data.Index65);
cada2f2dY = cada2f1dY;
cada2f2 = full(cada2f1);
cada2td1 = zeros(104438,1);
cada2td1(Gator2Data.Index78) = cada1td1dY;
cada2td1(Gator2Data.Index79) = cada2td1(Gator2Data.Index79) + cada2f2dY;
cada1td1dY = cada2td1;
cada1td1 = cada1td1 + cada2f2;
cada1f7dYdY = cada1td1dY; cada1f7dY = cada1td1;
cada1f7 = cada1f6*F.f;
cada1td1 =  zeros(26632,1);
cada1td1(Gator1Data.Index66) = cada1f5dY;
cada2f1 = cada1td1(Gator1Data.Index67);
cada2f2dY = cada1f7dYdY;
cada2f2 = cada2f1 + cada1f7dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index67) = cada2f2;
X_quad.dYdY = cada1td1dY; X_quad.dY = cada1td1;
X_quad.f = cada1f5 + cada1f7;
%User Line: X_quad=data.repXend_mat*X(1:2:end-1,:)+delta_t*data.AxHS*F;
cada1td1 = sparse(Gator1Data.Index68,Gator1Data.Index69,U.dY,199,199);
cada1td1 = data.AuHS*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index70);
U_quad.dY = full(cada2f1);
U_quad.f = data.AuHS*U.f;
%User Line: U_quad=data.AuHS*U;
cada1td1dY = F.dYdY(Gator2Data.Index80);
cada1td1 = sparse(Gator1Data.Index71,Gator1Data.Index72,F.dY,297,1098);
cada2td1 = sparse(Gator2Data.Index81,Gator2Data.Index82,cada1td1dY,297,3489);
cada2td1 = data.AfHS*cada2td1;
cada2td1 = cada2td1(:);
cada1td1dY = full(cada2td1(Gator2Data.Index83));
cada1td1 = data.AfHS*cada1td1;
cada1td1 = cada1td1(:);
cada2f1dY = cada1td1dY(Gator2Data.Index84);
cada2f1 = cada1td1(Gator1Data.Index73);
dX_quad.dYdY = cada2f1dY;
dX_quad.dY = full(cada2f1);
dX_quad.f = data.AfHS*F.f;
%User Line: dX_quad=data.AfHS*F;
P_quad.f = repmat(p,1386,1);
%User Line: P_quad=repmat(p,data.M_quad,1);
cada1tempdY = delta_t.dY(Gator1Data.Index74);
cada1tf1 = data.tau_quad(Gator1Data.Index76);
cada2f1 = cada1tf1(:);
cada2f2 = cada1tempdY(Gator1Data.Index75);
T_quad.dY = cada2f1.*cada2f2;
T_quad.f = data.tau_quad*delta_t.f;
%User Line: T_quad=data.tau_quad*delta_t;
cada1f1dYdY = X_quad.dYdY(Gator2Data.Index85);
cada1f1dY = X_quad.dY(Gator1Data.Index77);
cada1f1 = X_quad.f(:,1);
cada1f2dYdY = X_quad.dYdY(Gator2Data.Index86);
cada1f2dY = X_quad.dY(Gator1Data.Index78);
cada1f2 = X_quad.f(:,2);
cada1f3dYdY = X_quad.dYdY(Gator2Data.Index87);
cada1f3dY = X_quad.dY(Gator1Data.Index79);
cada1f3 = X_quad.f(:,3);
cada1f4dY = U_quad.dY(Gator1Data.Index80);
cada1f4 = U_quad.f(:,1);
cada1tf2dY = cada1f3dY(Gator2Data.Index88);
cada1tf2 = cada1f3(Gator1Data.Index81);
cada2f1dY = cada1tf2dY;
cada2f1 = cada1tf2(:);
cada2tf2 = cada2f1(Gator2Data.Index89);
cada2f2dY = 2.*cada2tf2(:).^(2-1).*cada2f1dY;
cada2f2 = cada2f1.^2;
cada2tf2 = cada2f2(Gator2Data.Index90);
cada2f3dY = -(-1)./cada2tf2(:).^2.*cada2f2dY;
cada2f3 = (-1)./cada2f2;
cada2tf1 = cada1f3dY(Gator2Data.Index91);
cada2td1 = cada2tf1(:).*cada2f3dY;
cada2tf1 = cada2f3(Gator2Data.Index92);
cada2td1(Gator2Data.Index93) = cada2td1(Gator2Data.Index93) + cada2tf1(:).*cada1f3dYdY;
cada1f5dYdY = cada2td1;
cada1f5dY = cada2f3.*cada1f3dY;
cada1f5 = 1./cada1f3;
cada1tf2dY = cada1f2dY(Gator2Data.Index94);
cada1tf2 = cada1f2(Gator1Data.Index82);
cada2f1dY = cada1tf2dY;
cada2f1 = cada1tf2(:);
cada2tf2 = cada2f1(Gator2Data.Index95);
cada2f2dY = 1.*cada2tf2(:).^(1-1).*cada2f1dY;
cada2f2 = cada2f1.^1;
cada2f3dY = 2.*cada2f2dY;
cada2f3 = 2*cada2f2;
cada2tf1 = cada1f2dY(Gator2Data.Index96);
cada2td1 = cada2tf1(:).*cada2f3dY;
cada2tf1 = cada2f3(Gator2Data.Index97);
cada2td1(Gator2Data.Index98) = cada2td1(Gator2Data.Index98) + cada2tf1(:).*cada1f2dYdY;
cada1f6dYdY = cada2td1;
cada1f6dY = cada2f3.*cada1f2dY;
cada1f6 = cada1f2.^2;
cada1f7dYdY = dyn_data.D0.*cada1f6dYdY;
cada1f7dY = dyn_data.D0*cada1f6dY;
cada1f7 = dyn_data.D0*cada1f6;
cada1f8dYdY = -cada1f1dYdY;
cada1f8dY = uminus(cada1f1dY);
cada1f8 = uminus(cada1f1);
cada1f9dYdY = cada1f8dYdY./dyn_data.H;
cada1f9dY = cada1f8dY/dyn_data.H;
cada1f9 = cada1f8/dyn_data.H;
cada1tf1dY = cada1f9dY(Gator2Data.Index99);
cada1tf1 = cada1f9(Gator1Data.Index83);
cada2f1dY = cada1tf1dY;
cada2f1 = cada1tf1(:);
cada2tf1 = cada2f1(Gator2Data.Index100);
cada2f2dY = exp(cada2tf1(:)).*cada2f1dY;
cada2f2 = exp(cada2f1);
cada2tf1 = cada1f9dY(Gator2Data.Index101);
cada2td1 = cada2tf1(:).*cada2f2dY;
cada2tf1 = cada2f2(Gator2Data.Index102);
cada2td1(Gator2Data.Index103) = cada2td1(Gator2Data.Index103) + cada2tf1(:).*cada1f9dYdY;
cada1f10dYdY = cada2td1;
cada1f10dY = cada2f2.*cada1f9dY;
cada1f10 = exp(cada1f9);
cada1tf1dY = cada1f10dY(Gator2Data.Index104);
cada1tf1 = cada1f10(Gator1Data.Index84);
cada1td1 =  zeros(13348,1);
cada2f1dY = cada1tf1dY;
cada2f1 = cada1tf1(:);
cada2tf1 = cada1f7dY(Gator2Data.Index105);
cada2td1 = zeros(105888,1);
cada2td1(Gator2Data.Index106) = cada2tf1(:).*cada2f1dY;
cada2tf1 = cada2f1(Gator2Data.Index107);
cada2td1(Gator2Data.Index108) = cada2td1(Gator2Data.Index108) + cada2tf1(:).*cada1f7dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f7dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index85) = cada2f2;
cada1tf1dY = cada1f7dY(Gator2Data.Index109);
cada1tf1 = cada1f7(Gator1Data.Index86);
cada2f1dY = cada1td1dY(Gator2Data.Index110);
cada2f1 = cada1td1(Gator1Data.Index87);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2tf1 = cada1f10dY(Gator2Data.Index111);
cada2td1 = zeros(79448,1);
cada2td1(Gator2Data.Index112) = cada2tf1(:).*cada2f2dY;
cada2tf1 = cada2f2(Gator2Data.Index113);
cada2td1(Gator2Data.Index114) = cada2td1(Gator2Data.Index114) + cada2tf1(:).*cada1f10dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f10dY;
cada2td1 = zeros(79448,1);
cada2td1(Gator2Data.Index115) = cada2f1dY;
cada2td1 = cada2td1 + cada2f3dY;
cada2f4dY = cada2td1;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(132456,1);
cada2td1(Gator2Data.Index116) = cada2f4dY;
cada2td1(Gator2Data.Index117) = cada1td1dY(Gator2Data.Index118);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index87) = cada2f4;
cada1f11dYdY = cada1td1dY; cada1f11dY = cada1td1;
cada1f11 = cada1f7.*cada1f10;
cada1td1 =  zeros(15986,1);
cada1td1(Gator1Data.Index88) = cada1f4dY;
cada2f1 = cada1td1(Gator1Data.Index89);
cada2f2dY = -cada1f11dYdY;
cada2f2 = uminus(cada1f11dY);
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada1td1dY = cada2f3dY;
cada1td1(Gator1Data.Index89) = cada2f3;
cada1f12dYdY = cada1td1dY; cada1f12dY = cada1td1;
cada1f12 = cada1f4 - cada1f11;
cada1tf1dY = cada1f12dY(Gator2Data.Index119);
cada1tf1 = cada1f12(Gator1Data.Index90);
cada1td1 =  zeros(18694,1);
cada2f1dY = cada1tf1dY;
cada2f1 = cada1tf1(:);
cada2tf1 = cada1f5dY(Gator2Data.Index120);
cada2td1 = zeros(110884,1);
cada2td1(Gator2Data.Index121) = cada2tf1(:).*cada2f1dY;
cada2tf1 = cada2f1(Gator2Data.Index122);
cada2td1(Gator2Data.Index123) = cada2td1(Gator2Data.Index123) + cada2tf1(:).*cada1f5dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f5dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index91) = cada2f2;
cada1tf1dY = cada1f5dY(Gator2Data.Index124);
cada1tf1 = cada1f5(Gator1Data.Index92);
cada2f1dY = cada1td1dY(Gator2Data.Index125);
cada2f1 = cada1td1(Gator1Data.Index93);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2tf1 = cada1f12dY(Gator2Data.Index126);
cada2td1 = zeros(174532,1);
cada2td1(Gator2Data.Index127) = cada2tf1(:).*cada2f2dY;
cada2tf1 = cada2f2(Gator2Data.Index128);
cada2td1(Gator2Data.Index129) = cada2td1(Gator2Data.Index129) + cada2tf1(:).*cada1f12dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f12dY;
cada2td1 = zeros(184828,1);
cada2td1(Gator2Data.Index130) = cada2f1dY;
cada2td1(Gator2Data.Index131) = cada2td1(Gator2Data.Index131) + cada2f3dY;
cada2f4dY = cada2td1;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(221960,1);
cada2td1(Gator2Data.Index132) = cada2f4dY;
cada2td1(Gator2Data.Index133) = cada1td1dY(Gator2Data.Index134);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index93) = cada2f4;
cada1f13dYdY = cada1td1dY; cada1f13dY = cada1td1;
cada1f13 = cada1f5.*cada1f12;
cada1f14dYdY = cada1f13dYdY; cada1f14dY = cada1f13dY;
cada1f14 = cada1f13 - dyn_data.grav;
cada1f15dY = uminus(cada1f4dY);
cada1f15 = uminus(cada1f4);
cada1f16dY = cada1f15dY/dyn_data.c;
cada1f16 = cada1f15/dyn_data.c;
cada1td1 =  zeros(33224,1);
cada1td1dY = cada1f2dYdY;
cada1td1(Gator1Data.Index94) = cada1f2dY;
cada2td1 = zeros(273518,1);
cada2td1(Gator2Data.Index135) = cada1f14dYdY;
cada2td1(Gator2Data.Index136) = cada1td1dY(Gator2Data.Index137);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index95) = cada1f14dY;
cada1td1(Gator1Data.Index96) = cada1f16dY;
Fp.dYdY = cada1td1dY; Fp.dY = cada1td1;
Fp.f = [cada1f2 cada1f14 cada1f16];
%User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
cada1td1 =  zeros(42548,1);
cada1td1dY = dX_quad.dYdY;
cada1td1(Gator1Data.Index97) = dX_quad.dY;
cada2f1dY = cada1td1dY(Gator2Data.Index138);
cada2f1 = cada1td1(Gator1Data.Index98);
cada2f2dY = -Fp.dYdY;
cada2f2 = uminus(Fp.dY);
cada2td1 = zeros(278806,1);
cada2td1(Gator2Data.Index139) = cada2f1dY;
cada2td1(Gator2Data.Index140) = cada2td1(Gator2Data.Index140) + cada2f2dY;
cada2f3dY = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(305246,1);
cada2td1(Gator2Data.Index141) = cada2f3dY;
cada2td1(Gator2Data.Index142) = cada1td1dY(Gator2Data.Index143);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index98) = cada2f3;
cada1f1dYdY = cada1td1dY; cada1f1dY = cada1td1;
cada1f1 = dX_quad.f - Fp.f;
cada1tf2dY = cada1f1dY(Gator2Data.Index144);
cada1tf2 = cada1f1(Gator1Data.Index99);
cada2f1dY = cada1tf2dY;
cada2f1 = cada1tf2(:);
cada2tf2 = cada2f1(Gator2Data.Index145);
cada2f2dY = 1.*cada2tf2(:).^(1-1).*cada2f1dY;
cada2f2 = cada2f1.^1;
cada2f3dY = 2.*cada2f2dY;
cada2f3 = 2*cada2f2;
cada2tf1 = cada1f1dY(Gator2Data.Index146);
cada2td1 = cada2tf1(:).*cada2f3dY;
cada2tf1 = cada2f3(Gator2Data.Index147);
cada2td1(Gator2Data.Index148) = cada2td1(Gator2Data.Index148) + cada2tf1(:).*cada1f1dYdY;
Res.dYdY = cada2td1;
Res.dY = cada2f3.*cada1f1dY;
Res.f = cada1f1.^2;
%User Line: Res=(dX_quad-Fp).^2;
cada1tempdY = delta_t.dY(Gator1Data.Index100);
cada1tf1 = data.DT_seg_node_mat(Gator1Data.Index102);
cada2f1 = cada1tf1(:);
cada2f2 = cada1tempdY(Gator1Data.Index101);
cada1f1dY = cada2f1.*cada2f2;
cada1f1 = delta_t.f*data.DT_seg_node_mat;
cada1f2dY = cada1f1dY/2;
cada1f2 = cada1f1/2;
cada1td1 = sparse(Gator1Data.Index103,Gator1Data.Index104,cada1f2dY,99,198);
cada2f1 = data.sum_nps_quad.';
cada1td1 = cada2f1*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index105);
cada1f3dY = full(cada2f1);
cada1f3 = cada1f2*data.sum_nps_quad;
cada1td2 = sparse(Gator1Data.Index106,Gator1Data.Index107,cada1f3dY,1386,198);
cada2f1dY = Res.dY(Gator2Data.Index149);
cada2f1 = Res.f.';
cada2td1 = sparse(Gator2Data.Index150,Gator2Data.Index151,cada2f1dY,1386,1796);
cada2td1 = cada1td2.'*cada2td1;
cada2td1 = cada2td1(:);
cada1td2dY = full(cada2td1(Gator2Data.Index152));
cada1td2 = cada2f1*cada1td2;
cada1td1 =  zeros(3168,1);
cada2f1dY = cada1td2dY(Gator2Data.Index153);
cada2f1 = cada1td2(Gator1Data.Index108);
cada1td1dY = cada2f1dY;
cada1td1(Gator1Data.Index109) = cada2f1;
cada1td2dY = Res.dYdY(Gator2Data.Index154);
cada1td2 = sparse(Gator1Data.Index110,Gator1Data.Index111,Res.dY,1386,1796);
cada2td2 = sparse(Gator2Data.Index155,Gator2Data.Index156,cada1f3dY,1386,198);
cada2td2 = cada1td2.'*cada2td2;
cada2td1 = zeros(35640,1);
cada2td1(Gator2Data.Index158) = cada2td2(Gator2Data.Index157);
cada2td2 = sparse(Gator2Data.Index159,Gator2Data.Index160,cada1td2dY,1386,28976);
cada2td2 = cada1f3*cada2td2;
cada2td2 = cada2td2(:);
cada2td1 = cada2td1 + full(cada2td2(Gator2Data.Index161));
cada1td2dY = cada2td1;
cada1td2 = cada1f3*cada1td2;
cada1td2 = cada1td2(:);
cada2f1dY = cada1td2dY(Gator2Data.Index162);
cada2f1 = cada1td2(Gator1Data.Index112);
cada2f2dY = cada2f1dY;
cada2f2 = full(cada2f1);
cada2td1 = zeros(35640,1);
cada2td1(Gator2Data.Index163) = cada1td1dY;
cada2td1 = cada2td1 + cada2f2dY;
cada1td1dY = cada2td1;
cada1td1 = cada1td1 + cada2f2;
Res_int.dYdY = cada1td1dY; Res_int.dY = cada1td1;
Res_int.f = cada1f3*Res.f;
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
%User Line: Fp=f(X_quad,U_quad,P_quad,T_quad,dyn_data);
%User Line: Res=(dX_quad-Fp).^2;
%User Line: Res_int=deltat*data.DT_seg_node_mat./2*data.sum_nps_quad*Res;
%User Line: % Compuation of integrated residual for each dynamics equation
cada1f1dYdY = Res_int.dYdY; cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1dY = cada1f1dYdY(Gator2Data.Index164);
cada1td1 = sparse(Gator1Data.Index113,Gator1Data.Index114,cada1f1dY,297,798);
cada2td1 = sparse(Gator2Data.Index165,Gator2Data.Index166,cada1td1dY,297,15876);
cada2td1 = data.ResConstScaleMat*cada2td1;
cada2td1 = cada2td1(:);
cada1td1dY = full(cada2td1(Gator2Data.Index167));
cada1td1 = data.ResConstScaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1dY = cada1td1dY(Gator2Data.Index168);
cada2f1 = cada1td1(Gator1Data.Index115);
Res_int_Const.dYdY = cada2f1dY;
Res_int_Const.dY = full(cada2f1);
Res_int_Const.f = data.ResConstScaleMat*cada1f1;
%User Line: Res_int_Const=data.ResConstScaleMat*Res_int(:);
Res_int_Const.dYdY = Res_int_Const.dYdY; Res_int_Const.dY = Res_int_Const.dY;
Res_int_Const.f = reshape(Res_int_Const.f,data.nps,n.f);
%User Line: Res_int_Const=reshape(Res_int_Const,data.nps,n);
cada2f1dY = Res_int_Const.dYdY(Gator2Data.Index169);
cada2f1 = sparse(Gator1Data.Index116,Gator1Data.Index117,Res_int_Const.dY,99,1796);
cada2td1 = sum(sparse(Gator2Data.Index170,Gator2Data.Index171,cada2f1dY,99,28976),1);
cada1td1dY = full(cada2td1(:));
cada1td1 = sum(cada2f1,1);
cada2f1dY = cada1td1dY;
cada2f1 = cada1td1(:);
cada1f1dYdY = cada2f1dY;
cada1f1dY = full(cada2f1);
cada1f1 = sum(Res_int_Const.f);
Res_intsum.dYdY = cada1f1dYdY; Res_intsum.dY = cada1f1dY;
Res_intsum.f = cada1f1.';
%User Line: Res_intsum=sum(Res_int_Const)';
%User Line: % Compuation of integrated residual norm
cada1f1dYdY = Res_int.dYdY; cada1f1dY = Res_int.dY;
cada1f1 = Res_int.f(:);
cada1td1dY = cada1f1dYdY(Gator2Data.Index172);
cada1td1 = sparse(Gator1Data.Index118,Gator1Data.Index119,cada1f1dY,297,798);
cada2td1 = sparse(Gator2Data.Index173,Gator2Data.Index174,cada1td1dY,297,15876);
cada2td1 = data.ResNormScaleMat*cada2td1;
cada2td1 = cada2td1(:);
cada1td1dY = full(cada2td1(Gator2Data.Index175));
cada1td1 = data.ResNormScaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1dY = cada1td1dY(Gator2Data.Index176);
cada2f1 = cada1td1(Gator1Data.Index120);
Res_int_Norm.dYdY = cada2f1dY;
Res_int_Norm.dY = full(cada2f1);
Res_int_Norm.f = data.ResNormScaleMat*cada1f1;
%User Line: Res_int_Norm=data.ResNormScaleMat*Res_int(:);
Res_int_Norm.dYdY = Res_int_Norm.dYdY; Res_int_Norm.dY = Res_int_Norm.dY;
Res_int_Norm.f = reshape(Res_int_Norm.f,data.nps,n.f);
%User Line: Res_int_Norm=reshape(Res_int_Norm,data.nps,n);
cada2f1dY = Res_int_Norm.dYdY(Gator2Data.Index177);
cada2f1 = sparse(Gator1Data.Index121,Gator1Data.Index122,Res_int_Norm.dY,99,1796);
cada2td1 = sum(sparse(Gator2Data.Index178,Gator2Data.Index179,cada2f1dY,99,28976),1);
cada1td1dY = full(cada2td1(:));
cada1td1 = sum(cada2f1,1);
cada2f1dY = cada1td1dY;
cada2f1 = cada1td1(:);
cada1f1dYdY = cada2f1dY;
cada1f1dY = full(cada2f1);
cada1f1 = sum(Res_int_Norm.f);
cada1td1 =  zeros(3,798);
cada1td1dY = cada1f1dYdY;
cada1td1(Gator1Data.Index123) = cada1f1dY;
cada2td1 = sum(sparse(Gator2Data.Index180,Gator2Data.Index181,cada1td1dY,3,15876),1);
cada1td1dY = full(cada2td1(:));
cada1td1 = sum(cada1td1,1);
ResNorm_intsum.dYdY = cada1td1dY;
ResNorm_intsum.dY = cada1td1(:);
ResNorm_intsum.f = sum(cada1f1);
%User Line: ResNorm_intsum=sum(sum(Res_int_Norm));
ResNorm_intsum.dY_size = 798;
ResNorm_intsum.dY_location = Gator1Data.Index124;
Res_intsum.dY_size = [3 798];
Res_intsum.dY_location = Gator1Data.Index125;
ResNorm_intsum.dYdY_size = [ResNorm_intsum.dY_size,798];
ResNorm_intsum.dYdY_location = [ResNorm_intsum.dY_location(Gator2Data.Index182,:), Gator2Data.Index183];
Res_intsum.dYdY_size = [Res_intsum.dY_size,798];
Res_intsum.dYdY_location = [Res_intsum.dY_location(Gator2Data.Index184,:), Gator2Data.Index185];
end


function ADiGator_LoadData()
global ADiGator_minresCost_YY
ADiGator_minresCost_YY = load('minresCost_YY.mat');
return
end