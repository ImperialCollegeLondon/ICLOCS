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

function dx = dynamics_YY(x,u,p,t,data)
global ADiGator_dynamics_YY
if isempty(ADiGator_dynamics_YY); ADiGator_LoadData(); end
Gator1Data = ADiGator_dynamics_YY.dynamics_YY.Gator1Data;
Gator2Data = ADiGator_dynamics_YY.dynamics_YY.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: % Vehicle Dynamics
%User Line: %
%User Line: % Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the BSD License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.0
%User Line: % 1 May 2018
%User Line: % iclocs@imperial.ac.uk
auxdata = data.auxdata;
%User Line: auxdata = data.auxdata;
v.dY = x.dY(:,3);
v.f = x.f(:,3);
%User Line: v = x(:,3);
theta.dY = x.dY(:,4);
theta.f = x.f(:,4);
%User Line: theta = x(:,4);
phi.dY = x.dY(:,5);
phi.f = x.f(:,5);
%User Line: phi = x(:,5);
a.dY = u.dY(:,1);
a.f = u.f(:,1);
%User Line: a=u(:,1);
u2.dY = u.dY(:,2);
u2.f = u.f(:,2);
%User Line: u2=u(:,2);
cada2f1dY = cos(theta.f).*theta.dY;
cada2f1 = sin(theta.f);
cada2f2dY = -cada2f1dY;
cada2f2 = uminus(cada2f1);
cada1f1dYdY = theta.dY.*cada2f2dY;
cada1f1dY = cada2f2.*theta.dY;
cada1f1 = cos(theta.f);
cada2f1 = size(v.dY,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dY = v.dY.*cada1f1dY;
cada2f1 = cada1f1.*v.dY;
cada1td1dY = cada2f1dY;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dY,1),2);
cada2td1(:,1) = cada1f1dY.*v.dY;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dYdY;
cada2f2dY = cada2td1;
cada2f2 = v.f.*cada1f1dY;
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index1) = cada2f3dY;
cada2td1(:,2) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada1td1(:,2) = cada2f3;
posx_dot.dYdY = cada1td1dY; posx_dot.dY = cada1td1;
posx_dot.f = v.f.*cada1f1;
%User Line: posx_dot=v.*cos(theta);
cada2f1dY = -sin(theta.f).*theta.dY;
cada2f1 = cos(theta.f);
cada1f1dYdY = theta.dY.*cada2f1dY;
cada1f1dY = cada2f1.*theta.dY;
cada1f1 = sin(theta.f);
cada2f1 = size(v.dY,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dY = v.dY.*cada1f1dY;
cada2f1 = cada1f1.*v.dY;
cada1td1dY = cada2f1dY;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dY,1),2);
cada2td1(:,1) = cada1f1dY.*v.dY;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dYdY;
cada2f2dY = cada2td1;
cada2f2 = v.f.*cada1f1dY;
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index2) = cada2f3dY;
cada2td1(:,2) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada1td1(:,2) = cada2f3;
posy_dot.dYdY = cada1td1dY; posy_dot.dY = cada1td1;
posy_dot.f = v.f.*cada1f1;
%User Line: posy_dot=v.*sin(theta);
v_dot.dY = a.dY;
v_dot.f = a.f;
%User Line: v_dot=a;
cada2f1dY = sec(phi.f).*tan(phi.f).*phi.dY;
cada2f1 = sec(phi.f);
cada2f2dY = 2.*cada2f1.^(2-1).*cada2f1dY;
cada2f2 = cada2f1.^2;
cada1f1dYdY = phi.dY.*cada2f2dY;
cada1f1dY = cada2f2.*phi.dY;
cada1f1 = tan(phi.f);
cada2f1 = size(v.dY,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dY = v.dY.*cada1f1dY;
cada2f1 = cada1f1.*v.dY;
cada1td1dY = cada2f1dY;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dY,1),2);
cada2td1(:,1) = cada1f1dY.*v.dY;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dYdY;
cada2f2dY = cada2td1;
cada2f2 = v.f.*cada1f1dY;
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index3) = cada2f3dY;
cada2td1(:,2) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f2dYdY = cada1td1dY; cada1f2dY = cada1td1;
cada1f2 = v.f.*cada1f1;
theta_dot.dYdY = cada1f2dYdY./auxdata.l_axes;
theta_dot.dY = cada1f2dY/auxdata.l_axes;
theta_dot.f = cada1f2/auxdata.l_axes;
%User Line: theta_dot=v.*tan(phi)./auxdata.l_axes;
phi_dot.dY = u2.dY;
phi_dot.f = u2.f;
%User Line: phi_dot=u2;
cada2f1 = size(posx_dot.f,1);
cada1td1 = zeros(cada2f1,8);
cada1td1dY = posx_dot.dYdY;
cada1td1(:,Gator1Data.Index1) = posx_dot.dY;
cada2td1 = zeros(size(cada1td1,1),6);
cada2td1(:,Gator2Data.Index4) = posy_dot.dYdY;
cada2td1(:,Gator2Data.Index5) = cada1td1dY(:,Gator2Data.Index6);
cada1td1dY = cada2td1;
cada1td1(:,Gator1Data.Index2) = posy_dot.dY;
cada1td1(:,7) = v_dot.dY;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index7) = theta_dot.dYdY;
cada2td1(:,Gator2Data.Index8) = cada1td1dY(:,Gator2Data.Index9);
cada1td1dY = cada2td1;
cada1td1(:,Gator1Data.Index3) = theta_dot.dY;
cada1td1(:,8) = phi_dot.dY;
dx.dYdY = cada1td1dY; dx.dY = cada1td1;
dx.f = [posx_dot.f posy_dot.f v_dot.f theta_dot.f phi_dot.f];
%User Line: dx = [posx_dot, posy_dot, v_dot, theta_dot, phi_dot];
dx.dY_size = [5 7];
dx.dY_location = Gator1Data.Index4;
dx.dYdY_size = [dx.dY_size,7];
dx.dYdY_location = [dx.dY_location(Gator2Data.Index10,:), Gator2Data.Index11];
end


function ADiGator_LoadData()
global ADiGator_dynamics_YY
ADiGator_dynamics_YY = load('dynamics_YY.mat');
return
end