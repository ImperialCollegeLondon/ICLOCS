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
%User Line: % Fed-batch fermentor Dynamics
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
x4.dY = x.dY(:,4);
x4.f = x.f(:,4);
%User Line: x4 = x(:,4);
u1.dY = u.dY(:,1);
u1.f = u.f(:,1);
%User Line: u1 = u(:,1);
cada1f1dY = 0.006.*x1.dY;
cada1f1 = 0.006*x1.f;
cada1td1 = zeros(size(cada1f1dY,1),2);
cada1td1(:,1) = cada1f1dY;
cada1td1(:,2) = cada1td1(:,2) + x3.dY;
cada1f2dY = cada1td1;
cada1f2 = cada1f1 + x3.f;
cada1td1 = zeros(size(x3.dY,1),2);
cada1td1(:,2) = x3.dY./cada1f2;
cada1tf1 = x3.f(:,Gator1Data.Index1);
cada1tf2 = cada1f2(:,Gator1Data.Index2);
cada1td1 = cada1td1 + -cada1tf1./cada1tf2.^2.*cada1f2dY;
cada1f3dY = cada1td1;
cada1f3 = x3.f./cada1f2;
h1.dY = 0.11.*cada1f3dY;
h1.f = 0.11*cada1f3;
%User Line: h1 = 0.11*(x3./(0.006*x1+x3));
cada1f1dY = 10.*x3.dY;
cada1f1 = 10*x3.f;
cada1f2dY = cada1f1dY;
cada1f2 = 1 + cada1f1;
cada1td1 = cada1f2.*x3.dY;
cada1td1 = cada1td1 + x3.f.*cada1f2dY;
cada1f3dY = cada1td1;
cada1f3 = x3.f.*cada1f2;
cada1f4dY = cada1f3dY;
cada1f4 = 0.0001 + cada1f3;
cada1td1 = x3.dY./cada1f4;
cada1td1 = cada1td1 + -x3.f./cada1f4.^2.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = x3.f./cada1f4;
h2.dY = 0.0055.*cada1f5dY;
h2.f = 0.0055*cada1f5;
%User Line: h2 = 0.0055*(x3./(0.0001+x3.*(1+10*x3)));
cada1tf1 = x1.f(:,Gator1Data.Index3);
cada1td1 = cada1tf1.*h1.dY;
cada1td1(:,1) = cada1td1(:,1) + h1.f.*x1.dY;
cada1f1dY = cada1td1;
cada1f1 = h1.f.*x1.f;
cada1f2dY = x1.dY./500;
cada1f2 = x1.f/500;
cada1td1 = zeros(size(cada1f2dY,1),2);
cada1td1(:,1) = cada1f2dY./x4.f;
cada1td1(:,2) = cada1td1(:,2) + -cada1f2./x4.f.^2.*x4.dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f2./x4.f;
cada1td1 = zeros(size(u1.dY,1),3);
cada1td1(:,3) = cada1f3.*u1.dY;
cada1tf1 = u1.f(:,Gator1Data.Index4);
cada1td1(:,Gator1Data.Index5) = cada1td1(:,Gator1Data.Index5) + cada1tf1.*cada1f3dY;
cada1f4dY = cada1td1;
cada1f4 = u1.f.*cada1f3;
cada1td1 = zeros(size(cada1f1dY,1),4);
cada1td1(:,Gator1Data.Index6) = cada1f1dY;
cada1td1(:,Gator1Data.Index7) = cada1td1(:,Gator1Data.Index7) + -cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f1 - cada1f4;
cada1temp1 = Gator1Data.Data1;
dx.dY = cada1f5dY;
dx.f = cada1temp1;
dx.f(:,1) = cada1f5;
%User Line: dx(:,1) = (h1.*x1-u1.*(x1./500./x4));
cada1td1 = zeros(size(h2.dY,1),2);
cada1td1(:,2) = x1.f.*h2.dY;
cada1td1(:,1) = cada1td1(:,1) + h2.f.*x1.dY;
cada1f1dY = cada1td1;
cada1f1 = h2.f.*x1.f;
cada1f2dY = 0.01.*x2.dY;
cada1f2 = 0.01*x2.f;
cada1td1 = zeros(size(cada1f1dY,1),3);
cada1td1(:,Gator1Data.Index8) = cada1f1dY;
cada1td1(:,2) = cada1td1(:,2) + -cada1f2dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f1 - cada1f2;
cada1f4dY = x2.dY./500;
cada1f4 = x2.f/500;
cada1td1 = zeros(size(cada1f4dY,1),2);
cada1td1(:,1) = cada1f4dY./x4.f;
cada1td1(:,2) = cada1td1(:,2) + -cada1f4./x4.f.^2.*x4.dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f4./x4.f;
cada1td1 = zeros(size(u1.dY,1),3);
cada1td1(:,3) = cada1f5.*u1.dY;
cada1tf1 = u1.f(:,Gator1Data.Index9);
cada1td1(:,Gator1Data.Index10) = cada1td1(:,Gator1Data.Index10) + cada1tf1.*cada1f5dY;
cada1f6dY = cada1td1;
cada1f6 = u1.f.*cada1f5;
cada1td1 = zeros(size(cada1f3dY,1),5);
cada1td1(:,Gator1Data.Index11) = cada1f3dY;
cada1td1(:,Gator1Data.Index12) = cada1td1(:,Gator1Data.Index12) + -cada1f6dY;
cada1f7dY = cada1td1;
cada1f7 = cada1f3 - cada1f6;
cada1td1 = zeros(size(dx.f,1),9);
cada1td1(:,Gator1Data.Index13) = cada1f7dY;
cada1td1(:,Gator1Data.Index14) = dx.dY(:,Gator1Data.Index15);
dx.dY = cada1td1;
dx.f(:,2) = cada1f7;
%User Line: dx(:,2) = (h2.*x1-0.01*x2-u1.*(x2./500./x4));
cada1f1dY = -h1.dY;
cada1f1 = uminus(h1.f);
cada1tf1 = x1.f(:,Gator1Data.Index16);
cada1td1 = cada1tf1.*cada1f1dY;
cada1td1(:,1) = cada1td1(:,1) + cada1f1.*x1.dY;
cada1f2dY = cada1td1;
cada1f2 = cada1f1.*x1.f;
cada1f3dY = cada1f2dY./0.47;
cada1f3 = cada1f2/0.47;
cada1td1 = zeros(size(h2.dY,1),2);
cada1td1(:,2) = x1.f.*h2.dY;
cada1td1(:,1) = cada1td1(:,1) + h2.f.*x1.dY;
cada1f4dY = cada1td1;
cada1f4 = h2.f.*x1.f;
cada1f5dY = cada1f4dY./1.2;
cada1f5 = cada1f4/1.2;
cada1td1 = cada1f3dY;
cada1td1 = cada1td1 + -cada1f5dY;
cada1f6dY = cada1td1;
cada1f6 = cada1f3 - cada1f5;
cada1f7dY = 0.029.*x3.dY;
cada1f7 = 0.029*x3.f;
cada1f8dY = x3.dY;
cada1f8 = 0.0001 + x3.f;
cada1td1 = cada1f7dY./cada1f8;
cada1td1 = cada1td1 + -cada1f7./cada1f8.^2.*cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = cada1f7./cada1f8;
cada1td1 = zeros(size(x1.dY,1),2);
cada1td1(:,1) = cada1f9.*x1.dY;
cada1td1(:,2) = cada1td1(:,2) + x1.f.*cada1f9dY;
cada1f10dY = cada1td1;
cada1f10 = x1.f.*cada1f9;
cada1td1 = cada1f6dY;
cada1td1 = cada1td1 + -cada1f10dY;
cada1f11dY = cada1td1;
cada1f11 = cada1f6 - cada1f10;
cada1td1 = zeros(size(u1.dY,1),2);
cada1td1(:,2) = u1.dY./x4.f;
cada1td1(:,1) = cada1td1(:,1) + -u1.f./x4.f.^2.*x4.dY;
cada1f12dY = cada1td1;
cada1f12 = u1.f./x4.f;
cada1f13dY = x3.dY./500;
cada1f13 = x3.f/500;
cada1f14dY = -cada1f13dY;
cada1f14 = 1 - cada1f13;
cada1tf1 = cada1f14(:,Gator1Data.Index17);
cada1td1 = zeros(size(cada1f12dY,1),3);
cada1td1(:,Gator1Data.Index18) = cada1tf1.*cada1f12dY;
cada1td1(:,1) = cada1td1(:,1) + cada1f12.*cada1f14dY;
cada1f15dY = cada1td1;
cada1f15 = cada1f12.*cada1f14;
cada1td1 = zeros(size(cada1f11dY,1),4);
cada1td1(:,Gator1Data.Index19) = cada1f11dY;
cada1td1(:,Gator1Data.Index20) = cada1td1(:,Gator1Data.Index20) + cada1f15dY;
cada1f16dY = cada1td1;
cada1f16 = cada1f11 + cada1f15;
cada1td1 = zeros(size(dx.f,1),13);
cada1td1(:,Gator1Data.Index21) = cada1f16dY;
cada1td1(:,Gator1Data.Index22) = dx.dY(:,Gator1Data.Index23);
dx.dY = cada1td1;
dx.f(:,3) = cada1f16;
%User Line: dx(:,3) = (-h1.*x1/0.47-h2.*x1/1.2-x1.*(0.029*x3./(0.0001+x3))+u1./x4.*(1-x3/500));
cada1f1dY = u1.dY./500;
cada1f1 = u1.f/500;
cada1td1 = zeros(size(dx.f,1),14);
cada1td1(:,14) = cada1f1dY;
cada1td1(:,Gator1Data.Index24) = dx.dY(:,Gator1Data.Index25);
dx.dY = cada1td1;
dx.f(:,4) = cada1f1;
%User Line: dx(:,4) = u1/500;
%User Line: %------------- END OF CODE --------------
dx.dY_size = [4,5];
dx.dY_location = Gator1Data.Index26;
end


function ADiGator_LoadData()
global ADiGator_dynamics_Y
ADiGator_dynamics_Y = load('dynamics_Y.mat');
return
end