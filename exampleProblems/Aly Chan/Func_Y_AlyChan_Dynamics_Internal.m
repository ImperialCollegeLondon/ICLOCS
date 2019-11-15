% This code was generated using ADiGator version 1.4
% ?2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function const = Func_Y_AlyChan_Dynamics_Internal(X_in,U,p,t,data)
global ADiGator_Func_Y_AlyChan_Dynamics_Internal
if isempty(ADiGator_Func_Y_AlyChan_Dynamics_Internal); ADiGator_LoadData(); end
Gator1Data = ADiGator_Func_Y_AlyChan_Dynamics_Internal.Func_Y_AlyChan_Dynamics_Internal.Gator1Data;
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
cada1td1 = zeros(397,1);
cada1td1(Gator1Data.Index6) = cada1f3dY;
cada1td1(Gator1Data.Index7) = cada1td1(Gator1Data.Index7) + cada1tempdY;
T.dY = cada1td1;
T.f = cada1f3 + t0.f;
%User Line: T=[0;data.tau_inc]*delta_t/ns+t0;
P.f = repmat(p,199,1);
%User Line: P=repmat(p,M,1);
cada1f1dY = X.dY(Gator1Data.Index8);
cada1f1 = X.f(:,1);
cada1f2dY = X.dY(Gator1Data.Index9);
cada1f2 = X.f(:,2);
cada1f3dY = U.dY(Gator1Data.Index10);
cada1f3 = U.f(:,1);
cada1temp1 = Gator1Data.Data1;
cada1f4dY = cada1f2dY;
cada1f4 = cada1temp1;
cada1f4(:,1) = cada1f2;
cada1td1 = zeros(398,1);
cada1td1(Gator1Data.Index11) = U.dY;
cada1td1(Gator1Data.Index12) = cada1f4dY(Gator1Data.Index13);
cada1f5dY = cada1td1;
cada1f5 = cada1f4;
cada1f5(:,2) = U.f;
cada1f6dY = 2.*cada1f2(:).^(2-1).*cada1f2dY;
cada1f6 = cada1f2.^2;
cada1f7dY = 2.*cada1f1(:).^(2-1).*cada1f1dY;
cada1f7 = cada1f1.^2;
cada1td1 = zeros(398,1);
cada1td1(Gator1Data.Index14) = cada1f6dY;
cada1td1(Gator1Data.Index15) = cada1td1(Gator1Data.Index15) + -cada1f7dY;
cada1f8dY = cada1td1;
cada1f8 = cada1f6 - cada1f7;
cada1f9dY = 0.5.*cada1f8dY;
cada1f9 = 0.5*cada1f8;
cada1td1 = zeros(796,1);
cada1td1(Gator1Data.Index16) = cada1f9dY;
cada1td1(Gator1Data.Index17) = cada1f5dY(Gator1Data.Index18);
dynF.dY = cada1td1;
dynF.f = cada1f5;
dynF.f(:,3) = cada1f9;
%User Line: dynF=f(X,U,P,T,vdat);
cada1f1dY = X.dY(Gator1Data.Index19);
cada1f1 = X.f.';
cada1f2 = n.f*M.f;
X_vect.dY = cada1f1dY;
X_vect.f = reshape(cada1f1,cada1f2,1);
%User Line: X_vect=reshape(X',n*M,1);
cada1f1dY = x0.dY;
cada1f1 = x0.f - data.x0t;
cada1f2dY = data.cx0.*cada1f1dY;
cada1f2 = cada1f1*data.cx0;
cada1td1 = sparse(Gator1Data.Index20,Gator1Data.Index21,X_vect.dY,597,597);
cada1td1 = data.map.A*cada1td1;
cada1td1 = cada1td1(:);
cada1f3dY = full(cada1td1(Gator1Data.Index22));
cada1f3 = data.map.A*X_vect.f;
cada1f4dY = dynF.dY(Gator1Data.Index23);
cada1f4 = dynF.f.';
cada1tempdY = delta_t.dY(Gator1Data.Index24);
cada1tf1 = cada1f4(Gator1Data.Index25);
cada1td1 = zeros(1990,1);
cada1td1(Gator1Data.Index26) = cada1tf1(:).*cada1tempdY;
cada1td1(Gator1Data.Index27) = cada1td1(Gator1Data.Index27) + delta_t.f.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = delta_t.f*cada1f4;
cada1f6 = M.f*n.f;
cada1f7dY = cada1f5dY;
cada1f7 = reshape(cada1f5,cada1f6,1);
cada1td1 = sparse(Gator1Data.Index28,Gator1Data.Index29,cada1f7dY,597,599);
cada1td1 = data.map.B*cada1td1;
cada1td1 = cada1td1(:);
cada1f8dY = full(cada1td1(Gator1Data.Index30));
cada1f8 = data.map.B*cada1f7;
cada1td1 = zeros(4653,1);
cada1td1(Gator1Data.Index31) = cada1f3dY;
cada1td1(Gator1Data.Index32) = cada1td1(Gator1Data.Index32) + cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = cada1f3 + cada1f8;
cada1td1 = zeros(4656,1);
cada1td1(Gator1Data.Index33) = cada1f2dY;
cada1td1(Gator1Data.Index34) = cada1f9dY;
const.dY = cada1td1;
const.f = [cada1f2;cada1f9];
%User Line: const=[(x0-data.x0t)*data.cx0;data.map.A*X_vect+data.map.B*reshape(delta_t*dynF',M*n,1)];
const.dY_size = [597,798];
const.dY_location = Gator1Data.Index35;
end


function ADiGator_LoadData()
global ADiGator_Func_Y_AlyChan_Dynamics_Internal
ADiGator_Func_Y_AlyChan_Dynamics_Internal = load('Func_Y_AlyChan_Dynamics_Internal.mat');
return
end