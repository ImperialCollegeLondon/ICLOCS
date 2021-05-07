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

function const = Func_YY_GoddardRocket_Dynamics_Internal(X_in,U,p,t,data)
global ADiGator_Func_YY_GoddardRocket_Dynamics_Internal
if isempty(ADiGator_Func_YY_GoddardRocket_Dynamics_Internal); ADiGator_LoadData(); end
Gator1Data = ADiGator_Func_YY_GoddardRocket_Dynamics_Internal.Func_YY_GoddardRocket_Dynamics_Internal.Gator1Data;
Gator2Data = ADiGator_Func_YY_GoddardRocket_Dynamics_Internal.Func_YY_GoddardRocket_Dynamics_Internal.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: %userFunction_Adigator_noConst_Scaling - Adigator template for user defined function (dynamics constraints) with direct collocation method (h-type) and automatic scaling
%User Line: %
%User Line: % Syntax:  [ const ] = userFunction_Adigator_noConst_Scaling( X_in,U,p,t,data)
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
X.dY = X_in.dY;
X.f = X_in.f;
%User Line: X=X_in;
cada1f1dY = X_in.dY(Gator1Data.Index1);
cada1f1 = X_in.f(1,:);
x0.dY = cada1f1dY;
x0.f = cada1f1.';
%User Line: x0=X_in(1,:)';
nt.f = data.sizes{1};
%User Line: nt=data.sizes{1};
n.f = data.sizes{3};
%User Line: n=data.sizes{3};
M.f = data.sizes{7};
%User Line: M=data.sizes{7};
ns.f = data.sizes{9};
%User Line: ns=data.sizes{9};
dimX1.f = size(X.f,1);
%User Line: dimX1=size(X,1);
dimX2.f = size(X.f,2);
%User Line: dimX2=size(X,2);
dimU1.f = size(U.f,1);
%User Line: dimU1=size(U,1);
dimU2.f = size(U.f,2);
%User Line: dimU2=size(U,2);
XshiftMat = data.scaling.XshiftMat;
%User Line: XshiftMat=data.scaling.XshiftMat;
UshiftMat = data.scaling.UshiftMat;
%User Line: UshiftMat=data.scaling.UshiftMat;
cada1f1dY = X.dY;
cada1f1 = X.f(:);
cada1td1 = sparse(Gator1Data.Index2,Gator1Data.Index3,cada1f1dY,597,597);
cada1td1 = data.scaling.XunscaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index4);
cada1f2dY = full(cada2f1);
cada1f2 = data.scaling.XunscaleMat*cada1f1;
cada1f3 = XshiftMat(:);
X.dY = cada1f2dY;
X.f = cada1f2 - cada1f3;
%User Line: X=data.scaling.XunscaleMat*X(:)-XshiftMat(:);
cada1f1dY = U.dY;
cada1f1 = U.f(:);
cada1td1 = sparse(Gator1Data.Index5,Gator1Data.Index6,cada1f1dY,199,199);
cada1td1 = data.scaling.UunscaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index7);
cada1f2dY = full(cada2f1);
cada1f2 = data.scaling.UunscaleMat*cada1f1;
cada1f3 = UshiftMat(:);
U.dY = cada1f2dY;
U.f = cada1f2 - cada1f3;
%User Line: U=data.scaling.UunscaleMat*U(:)-UshiftMat(:);
%User Line: % X=data.scaling.XunscaleMat*(X(:)-XshiftMat(:));
%User Line: % U=data.scaling.UunscaleMat*(U(:)-UshiftMat(:));
X.dY = X.dY;
X.f = reshape(X.f,dimX1.f,dimX2.f);
%User Line: X=reshape(X,dimX1,dimX2);
U.dY = U.dY;
U.f = reshape(U.f,dimU1.f,dimU2.f);
%User Line: U=reshape(U,dimU1,dimU2);
cadaconditional1 = isfield(data.data,'Pscale');
%User Line: cadaconditional1 = isfield(data.data,'Pscale');
%User Line: %     p=(p-data.data.Pshift)./data.data.Pscale;
%User Line: p=p./data.data.Pscale-data.data.Pshift;
t0.dY = t.dY(1);
t0.f = t.f(1);
%User Line: t0=t(1);
cada1f1 = length(t.f);
tf.dY = t.dY(2);
tf.f = t.f(cada1f1);
%User Line: tf=t(end);
cada1td1 =  zeros(2,1);
cada1td1(2) = tf.dY;
cada2f1 = cada1td1(1);
cada2f2 = uminus(t0.dY);
cada2f3 = cada2f1 + cada2f2;
cada1td1(1) = cada2f3;
delta_t.dY = cada1td1;
delta_t.f = tf.f - t0.f;
%User Line: delta_t=tf-t0;
cada1f1 = [0;data.tau_inc];
cada1tempdY = delta_t.dY(Gator1Data.Index8);
cada1tf1 = cada1f1(Gator1Data.Index10);
cada2f1 = cada1tf1(:);
cada2f2 = cada1tempdY(Gator1Data.Index9);
cada1f2dY = cada2f1.*cada2f2;
cada1f2 = cada1f1*delta_t.f;
cada1f3dY = cada1f2dY/ns.f;
cada1f3 = cada1f2/ns.f;
cada1tempdY = t0.dY(Gator1Data.Index11);
cada1td1 =  zeros(397,1);
cada1td1(Gator1Data.Index12) = cada1f3dY;
cada2f1 = cada1td1(Gator1Data.Index13);
cada2f2 = cada2f1 + cada1tempdY;
cada1td1(Gator1Data.Index13) = cada2f2;
T.dY = cada1td1;
T.f = cada1f3 + t0.f;
%User Line: T=[0;data.tau_inc]*delta_t/ns+t0;
P.f = repmat(p,199,1);
%User Line: P=repmat(p,M,1);
cada1f1dY = X.dY(Gator1Data.Index14);
cada1f1 = X.f(:,1);
cada1f2dY = X.dY(Gator1Data.Index15);
cada1f2 = X.f(:,2);
cada1f3dY = X.dY(Gator1Data.Index16);
cada1f3 = X.f(:,3);
cada1f4dY = U.dY(Gator1Data.Index17);
cada1f4 = U.f(:,1);
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
cada1f7dYdY = vdat.D0.*cada1f6dYdY;
cada1f7dY = vdat.D0*cada1f6dY;
cada1f7 = vdat.D0*cada1f6;
cada1f8dY = uminus(cada1f1dY);
cada1f8 = uminus(cada1f1);
cada1f9dY = cada1f8dY/vdat.H;
cada1f9 = cada1f8/vdat.H;
cada2f1dY = cada1f9dY;
cada2f1 = cada1f9(:);
cada2f2dY = exp(cada2f1(:)).*cada2f1dY;
cada2f2 = exp(cada2f1);
cada1f10dYdY = cada1f9dY(:).*cada2f2dY;
cada1f10dY = cada2f2.*cada1f9dY;
cada1f10 = exp(cada1f9);
cada1td1 =  zeros(398,1);
cada2f1dY = cada1f10dY;
cada2f1 = cada1f10(:);
cada2td1 = zeros(398,1);
cada2td1(Gator2Data.Index1) = cada1f7dY(:).*cada2f1dY;
cada2td1(Gator2Data.Index2) = cada2td1(Gator2Data.Index2) + cada2f1(:).*cada1f7dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f7dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index18) = cada2f2;
cada2f1 = cada1td1(Gator1Data.Index19);
cada2f2dY = cada1f7dY;
cada2f2 = cada1f7(:);
cada2td1 = zeros(398,1);
cada2td1(Gator2Data.Index3) = cada1f10dY(:).*cada2f2dY;
cada2td1(Gator2Data.Index4) = cada2td1(Gator2Data.Index4) + cada2f2(:).*cada1f10dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f10dY;
cada2f4dY = cada2f3dY;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(796,1);
cada2td1(Gator2Data.Index5) = cada2f4dY;
cada2td1(Gator2Data.Index6) = cada1td1dY(Gator2Data.Index7);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index19) = cada2f4;
cada1f11dYdY = cada1td1dY; cada1f11dY = cada1td1;
cada1f11 = cada1f7.*cada1f10;
cada1td1 =  zeros(597,1);
cada1td1(Gator1Data.Index20) = cada1f4dY;
cada2f1 = cada1td1(Gator1Data.Index21);
cada2f2dY = -cada1f11dYdY;
cada2f2 = uminus(cada1f11dY);
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada1td1dY = cada2f3dY;
cada1td1(Gator1Data.Index21) = cada2f3;
cada1f12dYdY = cada1td1dY; cada1f12dY = cada1td1;
cada1f12 = cada1f4 - cada1f11;
cada1td1 =  zeros(796,1);
cada2f1dY = cada1f12dY;
cada2f1 = cada1f12(:);
cada2tf1 = cada1f5dY(Gator2Data.Index8);
cada2td1 = zeros(796,1);
cada2td1(Gator2Data.Index9) = cada2tf1(:).*cada2f1dY;
cada2td1(Gator2Data.Index10) = cada2td1(Gator2Data.Index10) + cada2f1(:).*cada1f5dYdY;
cada2f2dY = cada2td1;
cada2f2 = cada2f1.*cada1f5dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index22) = cada2f2;
cada1tf1dY = cada1f5dY(Gator2Data.Index11);
cada1tf1 = cada1f5(Gator1Data.Index23);
cada2f1 = cada1td1(Gator1Data.Index24);
cada2f2dY = cada1tf1dY;
cada2f2 = cada1tf1(:);
cada2td1 = zeros(1393,1);
cada2td1(Gator2Data.Index12) = cada1f12dY(:).*cada2f2dY;
cada2tf1 = cada2f2(Gator2Data.Index13);
cada2td1(Gator2Data.Index14) = cada2td1(Gator2Data.Index14) + cada2tf1(:).*cada1f12dYdY;
cada2f3dY = cada2td1;
cada2f3 = cada2f2.*cada1f12dY;
cada2f4dY = cada2f3dY;
cada2f4 = cada2f1 + cada2f3;
cada2td1 = zeros(2189,1);
cada2td1(Gator2Data.Index15) = cada2f4dY;
cada2td1(Gator2Data.Index16) = cada1td1dY(Gator2Data.Index17);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index24) = cada2f4;
cada1f13dYdY = cada1td1dY; cada1f13dY = cada1td1;
cada1f13 = cada1f5.*cada1f12;
cada1f14dYdY = cada1f13dYdY; cada1f14dY = cada1f13dY;
cada1f14 = cada1f13 - vdat.grav;
cada1f15dY = uminus(cada1f4dY);
cada1f15 = uminus(cada1f4);
cada1f16dY = cada1f15dY/vdat.c;
cada1f16 = cada1f15/vdat.c;
cada1td1 =  zeros(1194,1);
cada1td1(Gator1Data.Index25) = cada1f2dY;
cada1td1dY = cada1f14dYdY;
cada1td1(Gator1Data.Index26) = cada1f14dY;
cada1td1(Gator1Data.Index27) = cada1f16dY;
dynF_org.dYdY = cada1td1dY; dynF_org.dY = cada1td1;
dynF_org.f = [cada1f2 cada1f14 cada1f16];
%User Line: [dynF_org]=f(X,U,P,T,vdat);
cada1f1dYdY = dynF_org.dYdY; cada1f1dY = dynF_org.dY;
cada1f1 = dynF_org.f(:);
cada1td1dY = cada1f1dYdY(Gator2Data.Index18);
cada1td1 = sparse(Gator1Data.Index28,Gator1Data.Index29,cada1f1dY,597,796);
cada2td1 = sparse(Gator2Data.Index19,Gator2Data.Index20,cada1td1dY,597,2189);
cada2td1 = data.scaling.XscaleMat*cada2td1;
cada2td1 = cada2td1(:);
cada1td1dY = full(cada2td1(Gator2Data.Index21));
cada1td1 = data.scaling.XscaleMat*cada1td1;
cada1td1 = cada1td1(:);
cada2f1dY = cada1td1dY(Gator2Data.Index22);
cada2f1 = cada1td1(Gator1Data.Index30);
dynF.dYdY = cada2f1dY;
dynF.dY = full(cada2f1);
dynF.f = data.scaling.XscaleMat*cada1f1;
%User Line: dynF=data.scaling.XscaleMat*dynF_org(:);
cada1f1dYdY = dynF.dYdY; cada1f1dY = dynF.dY;
cada1f1 = reshape(dynF.f,dimX1.f,dimX2.f);
dynF.dYdY = cada1f1dYdY(Gator2Data.Index23);
dynF.dY = cada1f1dY(Gator1Data.Index31);
dynF.f = cada1f1.';
%User Line: dynF=reshape(dynF,dimX1,dimX2)';
cada1f1dY = X_in.dY(Gator1Data.Index32);
cada1f1 = X_in.f.';
cada1f2 = n.f*M.f;
X_vect.dY = cada1f1dY;
X_vect.f = reshape(cada1f1,cada1f2,1);
%User Line: X_vect=reshape(X_in',n*M,1);
cada1f1dY = x0.dY;
cada1f1 = x0.f - data.x0t;
cada1f2dY = data.cx0*cada1f1dY;
cada1f2 = cada1f1*data.cx0;
cada1td1 = sparse(Gator1Data.Index33,Gator1Data.Index34,X_vect.dY,597,597);
cada1td1 = data.map.A*cada1td1;
cada1td1 = cada1td1(:);
cada2f1 = cada1td1(Gator1Data.Index35);
cada1f3dY = full(cada2f1);
cada1f3 = data.map.A*X_vect.f;
cada1tempdY = delta_t.dY(Gator1Data.Index36);
cada1tf1dY = dynF.dY(Gator2Data.Index24);
cada1tf1 = dynF.f(Gator1Data.Index37);
cada1td1 =  zeros(2388,1);
cada2f1dY = cada1tf1dY;
cada2f1 = cada1tf1(:);
cada2tf1 = cada1tempdY(Gator2Data.Index25);
cada2f2dY = cada2tf1(:).*cada2f1dY;
cada2f2 = cada2f1.*cada1tempdY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index38) = cada2f2;
cada2f1 = cada1td1(Gator1Data.Index39);
cada2tempdY = delta_t.dY(Gator2Data.Index26);
cada2tf1 = dynF.dY(Gator2Data.Index27);
cada2td1 = zeros(4577,1);
cada2td1(Gator2Data.Index28) = cada2tf1(:).*cada2tempdY;
cada2td1(Gator2Data.Index29) = cada2td1(Gator2Data.Index29) + delta_t.f.*dynF.dYdY;
cada2f2dY = cada2td1;
cada2f2 = delta_t.f*dynF.dY;
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(6965,1);
cada2td1(Gator2Data.Index30) = cada2f3dY;
cada2td1(Gator2Data.Index31) = cada1td1dY(Gator2Data.Index32);
cada1td1dY = cada2td1;
cada1td1(Gator1Data.Index39) = cada2f3;
cada1f4dYdY = cada1td1dY; cada1f4dY = cada1td1;
cada1f4 = delta_t.f*dynF.f;
cada1f5 = M.f*n.f;
cada1f6dYdY = cada1f4dYdY; cada1f6dY = cada1f4dY;
cada1f6 = reshape(cada1f4,cada1f5,1);
cada1td1dY = cada1f6dYdY(Gator2Data.Index33);
cada1td1 = sparse(Gator1Data.Index40,Gator1Data.Index41,cada1f6dY,597,798);
cada2td1 = sparse(Gator2Data.Index34,Gator2Data.Index35,cada1td1dY,597,5373);
cada2td1 = data.map.B*cada2td1;
cada2td1 = cada2td1(:);
cada1td1dY = full(cada2td1(Gator2Data.Index36));
cada1td1 = data.map.B*cada1td1;
cada1td1 = cada1td1(:);
cada2f1dY = cada1td1dY(Gator2Data.Index37);
cada2f1 = cada1td1(Gator1Data.Index42);
cada1f7dYdY = cada2f1dY;
cada1f7dY = full(cada2f1);
cada1f7 = data.map.B*cada1f6;
cada1td1 =  zeros(5247,1);
cada1td1(Gator1Data.Index43) = cada1f3dY;
cada2f1 = cada1td1(Gator1Data.Index44);
cada2f2dY = cada1f7dYdY;
cada2f2 = cada2f1 + cada1f7dY;
cada1td1dY = cada2f2dY;
cada1td1(Gator1Data.Index44) = cada2f2;
cada1f8dYdY = cada1td1dY; cada1f8dY = cada1td1;
cada1f8 = cada1f3 + cada1f7;
cada1td1 =  zeros(5250,1);
cada1td1(Gator1Data.Index45) = cada1f2dY;
cada1td1dY = cada1f8dYdY;
cada1td1(Gator1Data.Index46) = cada1f8dY;
const.dYdY = cada1td1dY; const.dY = cada1td1;
const.f = [cada1f2;cada1f8];
%User Line: const=[(x0-data.x0t)*data.cx0;data.map.A*X_vect+data.map.B*reshape(delta_t*dynF,M*n,1)];
const.dY_size = [597 798];
const.dY_location = Gator1Data.Index47;
const.dYdY_size = [const.dY_size,798];
const.dYdY_location = [const.dY_location(Gator2Data.Index38,:), Gator2Data.Index39];
end


function ADiGator_LoadData()
global ADiGator_Func_YY_GoddardRocket_Dynamics_Internal
ADiGator_Func_YY_GoddardRocket_Dynamics_Internal = load('Func_YY_GoddardRocket_Dynamics_Internal.mat');
return
end