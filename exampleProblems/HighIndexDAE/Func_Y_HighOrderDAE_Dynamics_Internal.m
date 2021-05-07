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

function const = Func_Y_HighOrderDAE_Dynamics_Internal(X_in,U,p,t,data)
global ADiGator_Func_Y_HighOrderDAE_Dynamics_Internal
if isempty(ADiGator_Func_Y_HighOrderDAE_Dynamics_Internal); ADiGator_LoadData(); end
Gator1Data = ADiGator_Func_Y_HighOrderDAE_Dynamics_Internal.Func_Y_HighOrderDAE_Dynamics_Internal.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %userFunction_Adigator_LGR_oneConst - Adigator template for user defined function (dynamics + equality or inequality path constraints) with direct collocation method (hp-type)
%User Line: %
%User Line: % Syntax:  [ const ] = userFunction_Adigator_LGR_oneConst( X_in,U,p,t,data)
%User Line: %
%User Line: % Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the MIT License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.5
%User Line: % 1 Aug 2019
%User Line: % iclocs@imperial.ac.uk
%User Line: %------------- BEGIN CODE --------------
f=data.data.InternalDynamics;
%User Line: f=data.data.InternalDynamics;
vdat=data.data;
%User Line: vdat=data.data;
X_Np1.dY = X_in.dY; X_Np1.f = X_in.f;
%User Line: X_Np1=X_in;
cada1f1 = size(X_Np1.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
X.dY = X_Np1.dY(Gator1Data.Index1);
X.f = X_Np1.f(cada1f3,:);
%User Line: X=X_Np1(1:end-1,:);
n.f = data.sizes{3};
%User Line: n=data.sizes{3};
ng.f = data.sizes{5};
%User Line: ng=data.sizes{5};
M.f = data.sizes{7};
%User Line: M=data.sizes{7};
D_structure = data.map.D_structure;
%User Line: D_structure=data.map.D_structure;
t_0.dY = t.dY(1);
t_0.f = t.f(1);
%User Line: t_0=t(1);
cada1f1 = length(t.f);
t_f.dY = t.dY(2);
t_f.f = t.f(cada1f1);
%User Line: t_f=t(end);
cada1td1 = zeros(2,1);
cada1td1(2) = t_f.dY;
cada1td1(1) = cada1td1(1) + -t_0.dY;
delta_t.dY = cada1td1;
delta_t.f = t_f.f - t_0.f;
%User Line: delta_t=t_f-t_0;
cada1f1dY = delta_t.dY./2;
cada1f1 = delta_t.f/2;
cada1tempdY = cada1f1dY(Gator1Data.Index2);
cada1tf1 = data.tau_inc(Gator1Data.Index4);
cada1f2dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index3);
cada1f2 = cada1f1*data.tau_inc;
cada1f3dY = delta_t.dY./2;
cada1f3 = delta_t.f/2;
cada1tempdY = cada1f3dY(Gator1Data.Index5);
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index6) = cada1f2dY;
cada1td1 = cada1td1 + cada1tempdY;
T.dY = cada1td1;
T.f = cada1f2 + cada1f3;
%User Line: T=delta_t/2*data.tau_inc+delta_t/2;
P.f = repmat(p,40,1);
%User Line: P=repmat(p,M,1);
cada1f1dY = X.dY(Gator1Data.Index7);
cada1f1 = X.f(:,1);
cada1f2dY = X.dY(Gator1Data.Index8);
cada1f2 = X.f(:,2);
cada1f3dY = X.dY(Gator1Data.Index9);
cada1f3 = X.f(:,3);
cada1f4dY = X.dY(Gator1Data.Index10);
cada1f4 = X.f(:,4);
cada1f5dY = U.dY(Gator1Data.Index11);
cada1f5 = U.f(:,1);
cada1f6dY = U.dY(Gator1Data.Index12);
cada1f6 = U.f(:,2);
cada1f7dY = -cada1f5dY;
cada1f7 = uminus(cada1f5);
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index13) = cada1f1(:).*cada1f7dY;
cada1td1(Gator1Data.Index14) = cada1td1(Gator1Data.Index14) + cada1f7(:).*cada1f1dY;
cada1f8dY = cada1td1;
cada1f8 = cada1f7.*cada1f1;
cada1f9dY = vdat.a.*cada1f2dY;
cada1f9 = vdat.a*cada1f2;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index15) = cada1f8dY;
cada1td1(Gator1Data.Index16) = cada1td1(Gator1Data.Index16) + -cada1f9dY;
cada1f10dY = cada1td1;
cada1f10 = cada1f8 - cada1f9;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index17) = cada1f3(:).*cada1f6dY;
cada1td1(Gator1Data.Index18) = cada1td1(Gator1Data.Index18) + cada1f6(:).*cada1f3dY;
cada1f11dY = cada1td1;
cada1f11 = cada1f6.*cada1f3;
cada1td1 = zeros(200,1);
cada1td1(Gator1Data.Index19) = cada1f10dY;
cada1td1(Gator1Data.Index20) = cada1td1(Gator1Data.Index20) + cada1f11dY;
cada1f12dY = cada1td1;
cada1f12 = cada1f10 + cada1f11;
cada1f13 = uminus(vdat.g);
cada1f14dY = 2.*cada1f5dY;
cada1f14 = 2*cada1f5;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index21) = cada1f3(:).*cada1f14dY;
cada1td1(Gator1Data.Index22) = cada1td1(Gator1Data.Index22) + cada1f14(:).*cada1f3dY;
cada1f15dY = cada1td1;
cada1f15 = cada1f14.*cada1f3;
cada1f16dY = -cada1f15dY;
cada1f16 = cada1f13 - cada1f15;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index23) = cada1f1(:).*cada1f6dY;
cada1td1(Gator1Data.Index24) = cada1td1(Gator1Data.Index24) + cada1f6(:).*cada1f1dY;
cada1f17dY = cada1td1;
cada1f17 = cada1f6.*cada1f1;
cada1td1 = zeros(160,1);
cada1td1(Gator1Data.Index25) = cada1f16dY;
cada1td1(Gator1Data.Index26) = cada1td1(Gator1Data.Index26) + -cada1f17dY;
cada1f18dY = cada1td1;
cada1f18 = cada1f16 - cada1f17;
cada1f19dY = vdat.a.*cada1f4dY;
cada1f19 = vdat.a*cada1f4;
cada1td1 = zeros(200,1);
cada1td1(Gator1Data.Index27) = cada1f18dY;
cada1td1(Gator1Data.Index28) = cada1td1(Gator1Data.Index28) + -cada1f19dY;
cada1f20dY = cada1td1;
cada1f20 = cada1f18 - cada1f19;
cada1f21dY = 2.*cada1f1(:).^(2-1).*cada1f1dY;
cada1f21 = cada1f1.^2;
cada1f22dY = 2.*cada1f3(:).^(2-1).*cada1f3dY;
cada1f22 = cada1f3.^2;
cada1td1 = zeros(80,1);
cada1td1(Gator1Data.Index29) = cada1f21dY;
cada1td1(Gator1Data.Index30) = cada1td1(Gator1Data.Index30) + cada1f22dY;
cada1f23dY = cada1td1;
cada1f23 = cada1f21 + cada1f22;
gConst.dY = cada1f23dY;
gConst.f = cada1f23 - vdat.L;
cada1td1 = zeros(480,1);
cada1td1(Gator1Data.Index31) = cada1f2dY;
cada1td1(Gator1Data.Index32) = cada1f12dY;
cada1td1(Gator1Data.Index33) = cada1f4dY;
cada1td1(Gator1Data.Index34) = cada1f20dY;
dynF.dY = cada1td1;
dynF.f = [cada1f2 cada1f12 cada1f4 cada1f20];
%User Line: [dynF,gConst]=f(X,U,P,T,vdat);
cada1f1 = M.f*ng.f;
g_all.dY = gConst.dY;
g_all.f = reshape(gConst.f,cada1f1,1);
%User Line: g_all=reshape(gConst,M*ng,1);
g_vect.dY = g_all.dY(Gator1Data.Index35);
g_vect.f = g_all.f(data.gAllidx);
%User Line: g_vect=g_all(data.gAllidx);
cada1td1 = zeros(2,1);
cada1td1(2) = t_f.dY;
cada1td1(1) = cada1td1(1) + -t_0.dY;
cada1f1dY = cada1td1;
cada1f1 = t_f.f - t_0.f;
cada1f2dY = cada1f1dY./2;
cada1f2 = cada1f1/2;
cada1f3 = data.tau_segment.';
cada1tempdY = cada1f2dY(Gator1Data.Index36);
cada1tf1 = cada1f3(Gator1Data.Index38);
cada1f4dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index37);
cada1f4 = cada1f2*cada1f3;
cada1td1 = zeros(2,1);
cada1td1(2) = t_f.dY;
cada1td1(1) = cada1td1(1) + t_0.dY;
cada1f5dY = cada1td1;
cada1f5 = t_f.f + t_0.f;
cada1f6dY = cada1f5dY./2;
cada1f6 = cada1f5/2;
cada1tempdY = cada1f6dY(Gator1Data.Index39);
cada1td1 = zeros(18,1);
cada1td1(Gator1Data.Index40) = cada1f4dY;
cada1td1 = cada1td1 + cada1tempdY;
t_segment.dY = cada1td1;
t_segment.f = cada1f4 + cada1f6;
%User Line: t_segment=(t_f-t_0)/2*data.tau_segment'+(t_f+t_0)/2;
cada1td1 = zeros(9,2);
cada1td1(Gator1Data.Index41) = t_segment.dY;
cada1td1 = data.t_segment_mat_m*cada1td1;
cada1td1 = cada1td1(:);
t_segment_end.dY = cada1td1(Gator1Data.Index42);
t_segment_end.f = data.t_segment_mat_m*t_segment.f;
%User Line: t_segment_end=data.t_segment_mat_m*t_segment;
cada1f1dY = t_segment_end.dY;
cada1f1 = diag(t_segment_end.f,0);
diag_t_segment_end.dY = cada1f1dY;
diag_t_segment_end.f = sparse(cada1f1);
%User Line: diag_t_segment_end=sparse(diag(t_segment_end));
cada1td1 = sparse(Gator1Data.Index43,Gator1Data.Index44,X_Np1.dY,41,164);
cada1td1 = D_structure*cada1td1;
cada1td1 = cada1td1(:);
cada1f1dY = full(cada1td1(Gator1Data.Index45));
cada1f1 = D_structure*X_Np1.f;
cada1td2 = sparse(Gator1Data.Index46,Gator1Data.Index47,diag_t_segment_end.dY,40,80);
cada1td2 = dynF.f.'*cada1td2;
cada1td1 = zeros(800,1);
cada1td1(Gator1Data.Index49) = cada1td2(Gator1Data.Index48);
cada1td2 = sparse(Gator1Data.Index50,Gator1Data.Index51,dynF.dY,40,480);
cada1td2 = diag_t_segment_end.f*cada1td2;
cada1td2 = cada1td2(:);
cada1td1(Gator1Data.Index53) = cada1td1(Gator1Data.Index53) + cada1td2(Gator1Data.Index52);
cada1f2dY = cada1td1;
cada1f2 = diag_t_segment_end.f*dynF.f;
cada1td1 = zeros(1680,1);
cada1td1(Gator1Data.Index54) = cada1f1dY;
cada1td1(Gator1Data.Index55) = cada1td1(Gator1Data.Index55) + -cada1f2dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f1 - cada1f2;
cada1f4 = M.f*n.f;
cada1f5dY = cada1f3dY;
cada1f5 = reshape(cada1f3,cada1f4,1);
cada1td1 = zeros(1760,1);
cada1td1(Gator1Data.Index56) = cada1f5dY;
cada1td1(Gator1Data.Index57) = g_vect.dY;
const.dY = cada1td1;
const.f = [cada1f5;g_vect.f];
%User Line: const=[reshape(D_structure*X_Np1-diag_t_segment_end*dynF,M*n,1);g_vect];
const.dY_size = [200,246];
const.dY_location = Gator1Data.Index58;
end


function ADiGator_LoadData()
global ADiGator_Func_Y_HighOrderDAE_Dynamics_Internal
ADiGator_Func_Y_HighOrderDAE_Dynamics_Internal = load('Func_Y_HighOrderDAE_Dynamics_Internal.mat');
return
end