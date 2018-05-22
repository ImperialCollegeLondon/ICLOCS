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

function dx = dynamics_Y(x,u,p,t,vdat)
global ADiGator_dynamics_Y
if isempty(ADiGator_dynamics_Y); ADiGator_LoadData(); end
Gator1Data = ADiGator_dynamics_Y.dynamics_Y.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Two-link Robot Arm Dynamics
%User Line: %
%User Line: % Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the BSD License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.0
%User Line: % 1 May 2018
%User Line: % iclocs@imperial.ac.uk
%User Line: %
%User Line: %------------- BEGIN CODE --------------
x1.dY = x.dY(:,1);
x1.f = x.f(:,1);
%User Line: x1 = x(:,1);
x2.dY = x.dY(:,2);
x2.f = x.f(:,2);
%User Line: x2 = x(:,2);
x3.dY = x.dY(:,3);
x3.f = x.f(:,3);
%User Line: x3 = x(:,3);
u1.dY = u.dY(:,1);
u1.f = u.f(:,1);
%User Line: u1 = u(:,1);
u2.dY = u.dY(:,2);
u2.f = u.f(:,2);
%User Line: u2 = u(:,2);
cada1f1dY = cos(x3.f).*x3.dY;
cada1f1 = sin(x3.f);
cada1f2dY = -sin(x3.f).*x3.dY;
cada1f2 = cos(x3.f);
cada1f3dY = 2.25.*cada1f2dY;
cada1f3 = 2.25*cada1f2;
cada1f4dY = 2.*x1.f.^(2-1).*x1.dY;
cada1f4 = x1.f.^2;
cada1td1 = zeros(size(cada1f3dY,1),2);
cada1td1(:,2) = cada1f4.*cada1f3dY;
cada1td1(:,1) = cada1td1(:,1) + cada1f3.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f3.*cada1f4;
cada1td1 = zeros(size(cada1f1dY,1),2);
cada1td1(:,2) = cada1f5.*cada1f1dY;
cada1tf1 = cada1f1(:,Gator1Data.Index1);
cada1td1 = cada1td1 + cada1tf1.*cada1f5dY;
cada1f6dY = cada1td1;
cada1f6 = cada1f1.*cada1f5;
cada1f7dY = 2.*x2.f.^(2-1).*x2.dY;
cada1f7 = x2.f.^2;
cada1f8dY = 2.*cada1f7dY;
cada1f8 = 2*cada1f7;
cada1td1 = zeros(size(cada1f6dY,1),3);
cada1td1(:,Gator1Data.Index2) = cada1f6dY;
cada1td1(:,2) = cada1td1(:,2) + cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = cada1f6 + cada1f8;
cada1td1 = zeros(size(u1.dY,1),2);
cada1td1(:,1) = u1.dY;
cada1td1(:,2) = cada1td1(:,2) + -u2.dY;
cada1f10dY = cada1td1;
cada1f10 = u1.f - u2.f;
cada1f11dY = 1.333333333333333.*cada1f10dY;
cada1f11 = 1.333333333333333*cada1f10;
cada1td1 = zeros(size(cada1f9dY,1),5);
cada1td1(:,Gator1Data.Index3) = cada1f9dY;
cada1td1(:,Gator1Data.Index4) = cada1td1(:,Gator1Data.Index4) + cada1f11dY;
cada1f12dY = cada1td1;
cada1f12 = cada1f9 + cada1f11;
cada1f13dY = -sin(x3.f).*x3.dY;
cada1f13 = cos(x3.f);
cada1f14dY = 1.5.*cada1f13dY;
cada1f14 = 1.5*cada1f13;
cada1td1 = zeros(size(cada1f14dY,1),2);
cada1td1(:,1) = u2.f.*cada1f14dY;
cada1td1(:,2) = cada1td1(:,2) + cada1f14.*u2.dY;
cada1f15dY = cada1td1;
cada1f15 = cada1f14.*u2.f;
cada1td1 = cada1f12dY;
cada1td1(:,Gator1Data.Index5) = cada1td1(:,Gator1Data.Index5) + -cada1f15dY;
cada1f16dY = cada1td1;
cada1f16 = cada1f12 - cada1f15;
cada1f17dY = cos(x3.f).*x3.dY;
cada1f17 = sin(x3.f);
cada1f18dY = 2.*cada1f17.^(2-1).*cada1f17dY;
cada1f18 = cada1f17.^2;
cada1f19dY = 2.25.*cada1f18dY;
cada1f19 = 2.25*cada1f18;
cada1f20dY = cada1f19dY;
cada1f20 = 0.8611111111111112 + cada1f19;
cada1tf1 = cada1f20(:,Gator1Data.Index6);
cada1td1 = cada1f16dY./cada1tf1;
cada1td1(:,3) = cada1td1(:,3) + -cada1f16./cada1f20.^2.*cada1f20dY;
cada1f21dY = cada1td1;
cada1f21 = cada1f16./cada1f20;
cada1temp1 = Gator1Data.Data1;
dx.dY = cada1f21dY;
dx.f = cada1temp1;
dx.f(:,1) = cada1f21;
%User Line: dx(:,1) = ( sin(x3).*(9/4*cos(x3).*x1.^2)+2*x2.^2 + 4/3*(u1-u2)           - 3/2*cos(x3).*u2 )./ (31/36 + 9/4*sin(x3).^2);
cada1f1dY = cos(x3.f).*x3.dY;
cada1f1 = sin(x3.f);
cada1f2dY = -sin(x3.f).*x3.dY;
cada1f2 = cos(x3.f);
cada1f3dY = 2.25.*cada1f2dY;
cada1f3 = 2.25*cada1f2;
cada1f4dY = 2.*x2.f.^(2-1).*x2.dY;
cada1f4 = x2.f.^2;
cada1td1 = zeros(size(cada1f3dY,1),2);
cada1td1(:,2) = cada1f4.*cada1f3dY;
cada1td1(:,1) = cada1td1(:,1) + cada1f3.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f3.*cada1f4;
cada1td1 = zeros(size(cada1f1dY,1),2);
cada1td1(:,2) = cada1f5.*cada1f1dY;
cada1tf1 = cada1f1(:,Gator1Data.Index7);
cada1td1 = cada1td1 + cada1tf1.*cada1f5dY;
cada1f6dY = cada1td1;
cada1f6 = cada1f1.*cada1f5;
cada1f7dY = 2.*x1.f.^(2-1).*x1.dY;
cada1f7 = x1.f.^2;
cada1f8dY = 3.5.*cada1f7dY;
cada1f8 = 3.5*cada1f7;
cada1td1 = zeros(size(cada1f6dY,1),3);
cada1td1(:,Gator1Data.Index8) = cada1f6dY;
cada1td1(:,1) = cada1td1(:,1) + cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = cada1f6 + cada1f8;
cada1f10dY = 2.333333333333334.*u2.dY;
cada1f10 = 2.333333333333334*u2.f;
cada1td1 = zeros(size(cada1f9dY,1),4);
cada1td1(:,Gator1Data.Index9) = cada1f9dY;
cada1td1(:,4) = cada1td1(:,4) + -cada1f10dY;
cada1f11dY = cada1td1;
cada1f11 = cada1f9 - cada1f10;
cada1f12dY = -sin(x3.f).*x3.dY;
cada1f12 = cos(x3.f);
cada1f13dY = 1.5.*cada1f12dY;
cada1f13 = 1.5*cada1f12;
cada1td1 = zeros(size(u1.dY,1),2);
cada1td1(:,1) = u1.dY;
cada1td1(:,2) = cada1td1(:,2) + -u2.dY;
cada1f14dY = cada1td1;
cada1f14 = u1.f - u2.f;
cada1td1 = zeros(size(cada1f13dY,1),3);
cada1td1(:,1) = cada1f14.*cada1f13dY;
cada1tf1 = cada1f13(:,Gator1Data.Index10);
cada1td1(:,Gator1Data.Index11) = cada1td1(:,Gator1Data.Index11) + cada1tf1.*cada1f14dY;
cada1f15dY = cada1td1;
cada1f15 = cada1f13.*cada1f14;
cada1td1 = zeros(size(cada1f11dY,1),5);
cada1td1(:,Gator1Data.Index12) = cada1f11dY;
cada1td1(:,Gator1Data.Index13) = cada1td1(:,Gator1Data.Index13) + cada1f15dY;
cada1f16dY = cada1td1;
cada1f16 = cada1f11 + cada1f15;
cada1f17dY = -cada1f16dY;
cada1f17 = uminus(cada1f16);
cada1f18dY = cos(x3.f).*x3.dY;
cada1f18 = sin(x3.f);
cada1f19dY = 2.*cada1f18.^(2-1).*cada1f18dY;
cada1f19 = cada1f18.^2;
cada1f20dY = 2.25.*cada1f19dY;
cada1f20 = 2.25*cada1f19;
cada1f21dY = cada1f20dY;
cada1f21 = 0.8611111111111112 + cada1f20;
cada1tf1 = cada1f21(:,Gator1Data.Index14);
cada1td1 = cada1f17dY./cada1tf1;
cada1td1(:,3) = cada1td1(:,3) + -cada1f17./cada1f21.^2.*cada1f21dY;
cada1f22dY = cada1td1;
cada1f22 = cada1f17./cada1f21;
cada1td1 = zeros(size(dx.f,1),10);
cada1td1(:,Gator1Data.Index15) = cada1f22dY;
cada1td1(:,Gator1Data.Index16) = dx.dY(:,Gator1Data.Index17);
dx.dY = cada1td1;
dx.f(:,2) = cada1f22;
%User Line: dx(:,2) = -( sin(x3).*(9/4*cos(x3).*x2.^2)+7/2*x1.^2 - 7/3*u2           + 3/2*cos(x3).*(u1-u2) )./ (31/36 + 9/4*sin(x3).^2);
cada1td1 = zeros(size(x2.dY,1),2);
cada1td1(:,2) = x2.dY;
cada1td1(:,1) = cada1td1(:,1) + -x1.dY;
cada1f1dY = cada1td1;
cada1f1 = x2.f - x1.f;
cada1td1 = zeros(size(dx.f,1),12);
cada1td1(:,Gator1Data.Index18) = cada1f1dY;
cada1td1(:,Gator1Data.Index19) = dx.dY(:,Gator1Data.Index20);
dx.dY = cada1td1;
dx.f(:,3) = cada1f1;
%User Line: dx(:,3) = x2-x1;
cada1td1 = zeros(size(dx.f,1),13);
cada1td1(:,4) = x1.dY;
cada1td1(:,Gator1Data.Index21) = dx.dY(:,Gator1Data.Index22);
dx.dY = cada1td1;
dx.f(:,4) = x1.f;
%User Line: dx(:,4) = x1;
%User Line: %------------- END OF CODE --------------
dx.dY_size = [4,6];
dx.dY_location = Gator1Data.Index23;
end


function ADiGator_LoadData()
global ADiGator_dynamics_Y
ADiGator_dynamics_Y = load('dynamics_Y.mat');
return
end