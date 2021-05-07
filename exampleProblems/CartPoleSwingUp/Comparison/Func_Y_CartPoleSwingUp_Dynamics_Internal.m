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

function const = Func_Y_CartPoleSwingUp_Dynamics_Internal(X_in,U,p,t,data)
global ADiGator_Func_Y_CartPoleSwingUp_Dynamics_Internal
if isempty(ADiGator_Func_Y_CartPoleSwingUp_Dynamics_Internal); ADiGator_LoadData(); end
Gator1Data = ADiGator_Func_Y_CartPoleSwingUp_Dynamics_Internal.Func_Y_CartPoleSwingUp_Dynamics_Internal.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %userFunction_Adigator_noConst - Adigator template for user defined function (dynamics constraints) with direct collocation method (h-type)
%User Line: %
%User Line: % Syntax:  [ const ] = userFunction_Adigator_noConst( X_in,U,p,t,data)
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
X.dY = X_in.dY; X.f = X_in.f;
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
t0.dY = t.dY(1);
t0.f = t.f(1);
%User Line: t0=t(1);
cada1f1 = length(t.f);
tf.dY = t.dY(2);
tf.f = t.f(cada1f1);
%User Line: tf=t(end);
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dY;
cada1td1(1) = cada1td1(1) + -t0.dY;
delta_t.dY = cada1td1;
delta_t.f = tf.f - t0.f;
%User Line: delta_t=tf-t0;
cada1f1 = [0;data.tau_inc];
cada1tempdY = delta_t.dY(Gator1Data.Index2);
cada1tf1 = cada1f1(Gator1Data.Index4);
cada1f2dY = cada1tf1(:).*cada1tempdY(Gator1Data.Index3);
cada1f2 = cada1f1*delta_t.f;
cada1f3dY = cada1f2dY./ns.f;
cada1f3 = cada1f2/ns.f;
cada1tempdY = t0.dY(Gator1Data.Index5);
cada1td1 = zeros(29,1);
cada1td1(Gator1Data.Index6) = cada1f3dY;
cada1td1(Gator1Data.Index7) = cada1td1(Gator1Data.Index7) + cada1tempdY;
T.dY = cada1td1;
T.f = cada1f3 + t0.f;
%User Line: T=[0;data.tau_inc]*delta_t/ns+t0;
P.f = repmat(p,15,1);
%User Line: P=repmat(p,M,1);
cada1f1dY = X.dY(Gator1Data.Index8);
cada1f1 = X.f(:,1);
cada1f2dY = X.dY(Gator1Data.Index9);
cada1f2 = X.f(:,2);
cada1f3dY = X.dY(Gator1Data.Index10);
cada1f3 = X.f(:,3);
cada1f4dY = X.dY(Gator1Data.Index11);
cada1f4 = X.f(:,4);
cada1f5dY = U.dY(Gator1Data.Index12);
cada1f5 = U.f(:,1);
cada1f6 = vdat.L*vdat.m2;
cada1f7dY = cos(cada1f4(:)).*cada1f4dY;
cada1f7 = sin(cada1f4);
cada1f8dY = cada1f6.*cada1f7dY;
cada1f8 = cada1f6*cada1f7;
cada1f9dY = 2.*cada1f2(:).^(2-1).*cada1f2dY;
cada1f9 = cada1f2.^2;
cada1td1 = zeros(30,1);
cada1td1(Gator1Data.Index13) = cada1f9(:).*cada1f8dY;
cada1td1(Gator1Data.Index14) = cada1td1(Gator1Data.Index14) + cada1f8(:).*cada1f9dY;
cada1f10dY = cada1td1;
cada1f10 = cada1f8.*cada1f9;
cada1td1 = zeros(45,1);
cada1td1(Gator1Data.Index15) = cada1f10dY;
cada1td1(Gator1Data.Index16) = cada1td1(Gator1Data.Index16) + cada1f5dY;
cada1f11dY = cada1td1;
cada1f11 = cada1f10 + cada1f5;
cada1f12 = vdat.m2*vdat.g;
cada1f13dY = -sin(cada1f4(:)).*cada1f4dY;
cada1f13 = cos(cada1f4);
cada1f14dY = cada1f12.*cada1f13dY;
cada1f14 = cada1f12*cada1f13;
cada1f15dY = cos(cada1f4(:)).*cada1f4dY;
cada1f15 = sin(cada1f4);
cada1td1 = cada1f15(:).*cada1f14dY;
cada1td1 = cada1td1 + cada1f14(:).*cada1f15dY;
cada1f16dY = cada1td1;
cada1f16 = cada1f14.*cada1f15;
cada1td1 = cada1f11dY;
cada1td1(Gator1Data.Index17) = cada1td1(Gator1Data.Index17) + cada1f16dY;
cada1f17dY = cada1td1;
cada1f17 = cada1f11 + cada1f16;
cada1f18dY = -sin(cada1f4(:)).*cada1f4dY;
cada1f18 = cos(cada1f4);
cada1f19dY = 2.*cada1f18(:).^(2-1).*cada1f18dY;
cada1f19 = cada1f18.^2;
cada1f20dY = -cada1f19dY;
cada1f20 = 1 - cada1f19;
cada1f21dY = vdat.m2.*cada1f20dY;
cada1f21 = vdat.m2*cada1f20;
cada1f22dY = cada1f21dY;
cada1f22 = vdat.m1 + cada1f21;
cada1tf1 = cada1f22(Gator1Data.Index18);
cada1td1 = cada1f17dY./cada1tf1(:);
cada1td1(Gator1Data.Index19) = cada1td1(Gator1Data.Index19) + -cada1f17(:)./cada1f22(:).^2.*cada1f22dY;
cada1f23dY = cada1td1;
cada1f23 = cada1f17./cada1f22;
cada1temp1 = Gator1Data.Data1;
cada1f24dY = cada1f23dY;
cada1f24 = cada1temp1;
cada1f24(:,1) = cada1f23;
cada1f25 = vdat.L*vdat.m2;
cada1f26dY = -sin(cada1f4(:)).*cada1f4dY;
cada1f26 = cos(cada1f4);
cada1f27dY = cada1f25.*cada1f26dY;
cada1f27 = cada1f25*cada1f26;
cada1f28dY = cos(cada1f4(:)).*cada1f4dY;
cada1f28 = sin(cada1f4);
cada1td1 = cada1f28(:).*cada1f27dY;
cada1td1 = cada1td1 + cada1f27(:).*cada1f28dY;
cada1f29dY = cada1td1;
cada1f29 = cada1f27.*cada1f28;
cada1f30dY = 2.*cada1f2(:).^(2-1).*cada1f2dY;
cada1f30 = cada1f2.^2;
cada1td1 = zeros(30,1);
cada1td1(Gator1Data.Index20) = cada1f30(:).*cada1f29dY;
cada1td1(Gator1Data.Index21) = cada1td1(Gator1Data.Index21) + cada1f29(:).*cada1f30dY;
cada1f31dY = cada1td1;
cada1f31 = cada1f29.*cada1f30;
cada1f32dY = -sin(cada1f4(:)).*cada1f4dY;
cada1f32 = cos(cada1f4);
cada1td1 = zeros(30,1);
cada1td1(Gator1Data.Index22) = cada1f32(:).*cada1f5dY;
cada1td1(Gator1Data.Index23) = cada1td1(Gator1Data.Index23) + cada1f5(:).*cada1f32dY;
cada1f33dY = cada1td1;
cada1f33 = cada1f5.*cada1f32;
cada1td1 = zeros(45,1);
cada1td1(Gator1Data.Index24) = cada1f31dY;
cada1td1(Gator1Data.Index25) = cada1td1(Gator1Data.Index25) + cada1f33dY;
cada1f34dY = cada1td1;
cada1f34 = cada1f31 + cada1f33;
cada1f35 = vdat.m1 + vdat.m2;
cada1f36 = cada1f35*vdat.g;
cada1f37dY = cos(cada1f4(:)).*cada1f4dY;
cada1f37 = sin(cada1f4);
cada1f38dY = cada1f36.*cada1f37dY;
cada1f38 = cada1f36*cada1f37;
cada1td1 = cada1f34dY;
cada1td1(Gator1Data.Index26) = cada1td1(Gator1Data.Index26) + cada1f38dY;
cada1f39dY = cada1td1;
cada1f39 = cada1f34 + cada1f38;
cada1f40dY = -cada1f39dY;
cada1f40 = uminus(cada1f39);
cada1f41 = vdat.L*vdat.m1;
cada1f42 = vdat.L*vdat.m2;
cada1f43dY = -sin(cada1f4(:)).*cada1f4dY;
cada1f43 = cos(cada1f4);
cada1f44dY = 2.*cada1f43(:).^(2-1).*cada1f43dY;
cada1f44 = cada1f43.^2;
cada1f45dY = -cada1f44dY;
cada1f45 = 1 - cada1f44;
cada1f46dY = cada1f42.*cada1f45dY;
cada1f46 = cada1f42*cada1f45;
cada1f47dY = cada1f46dY;
cada1f47 = cada1f41 + cada1f46;
cada1tf1 = cada1f47(Gator1Data.Index27);
cada1td1 = cada1f40dY./cada1tf1(:);
cada1td1(Gator1Data.Index28) = cada1td1(Gator1Data.Index28) + -cada1f40(:)./cada1f47(:).^2.*cada1f47dY;
cada1f48dY = cada1td1;
cada1f48 = cada1f40./cada1f47;
cada1td1 = zeros(90,1);
cada1td1(Gator1Data.Index29) = cada1f48dY;
cada1td1(Gator1Data.Index30) = cada1f24dY(Gator1Data.Index31);
cada1f49dY = cada1td1;
cada1f49 = cada1f24;
cada1f49(:,2) = cada1f48;
cada1td1 = zeros(105,1);
cada1td1(Gator1Data.Index32) = cada1f1dY;
cada1td1(Gator1Data.Index33) = cada1f49dY(Gator1Data.Index34);
cada1f50dY = cada1td1;
cada1f50 = cada1f49;
cada1f50(:,3) = cada1f1;
cada1td1 = zeros(120,1);
cada1td1(Gator1Data.Index35) = cada1f2dY;
cada1td1(Gator1Data.Index36) = cada1f50dY(Gator1Data.Index37);
dynF.dY = cada1td1;
dynF.f = cada1f50;
dynF.f(:,4) = cada1f2;
%User Line: dynF=f(X,U,P,T,vdat);
cada1f1dY = X.dY(Gator1Data.Index38);
cada1f1 = X.f.';
cada1f2 = n.f*M.f;
X_vect.dY = cada1f1dY;
X_vect.f = reshape(cada1f1,cada1f2,1);
%User Line: X_vect=reshape(X',n*M,1);
cada1f1dY = x0.dY;
cada1f1 = x0.f - data.x0t;
cada1f2dY = data.cx0.*cada1f1dY;
cada1f2 = cada1f1*data.cx0;
cada1td1 = sparse(Gator1Data.Index39,Gator1Data.Index40,X_vect.dY,60,60);
cada1td1 = data.map.A*cada1td1;
cada1td1 = cada1td1(:);
cada1f3dY = full(cada1td1(Gator1Data.Index41));
cada1f3 = data.map.A*X_vect.f;
cada1f4dY = dynF.dY(Gator1Data.Index42);
cada1f4 = dynF.f.';
cada1tempdY = delta_t.dY(Gator1Data.Index43);
cada1tf1 = cada1f4(Gator1Data.Index44);
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index45) = cada1tf1(:).*cada1tempdY;
cada1td1(Gator1Data.Index46) = cada1td1(Gator1Data.Index46) + delta_t.f.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = delta_t.f*cada1f4;
cada1f6 = M.f*n.f;
cada1f7dY = cada1f5dY;
cada1f7 = reshape(cada1f5,cada1f6,1);
cada1td1 = sparse(Gator1Data.Index47,Gator1Data.Index48,cada1f7dY,60,62);
cada1td1 = data.map.B*cada1td1;
cada1td1 = cada1td1(:);
cada1f8dY = full(cada1td1(Gator1Data.Index49));
cada1f8 = data.map.B*cada1f7;
cada1td1 = zeros(504,1);
cada1td1(Gator1Data.Index50) = cada1f3dY;
cada1td1(Gator1Data.Index51) = cada1td1(Gator1Data.Index51) + cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = cada1f3 + cada1f8;
cada1td1 = zeros(508,1);
cada1td1(Gator1Data.Index52) = cada1f2dY;
cada1td1(Gator1Data.Index53) = cada1f9dY;
const.dY = cada1td1;
const.f = [cada1f2;cada1f9];
%User Line: const=[(x0-data.x0t)*data.cx0;data.map.A*X_vect+data.map.B*reshape(delta_t*dynF',M*n,1)];
const.dY_size = [60,77];
const.dY_location = Gator1Data.Index54;
end


function ADiGator_LoadData()
global ADiGator_Func_Y_CartPoleSwingUp_Dynamics_Internal
ADiGator_Func_Y_CartPoleSwingUp_Dynamics_Internal = load('Func_Y_CartPoleSwingUp_Dynamics_Internal.mat');
return
end