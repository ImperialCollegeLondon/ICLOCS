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

function dx = dynamics_YY(x,u,p,t,vdat)
global ADiGator_dynamics_YY
if isempty(ADiGator_dynamics_YY); ADiGator_LoadData(); end
Gator1Data = ADiGator_dynamics_YY.dynamics_YY.Gator1Data;
Gator2Data = ADiGator_dynamics_YY.dynamics_YY.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: %Double Integrator Dynamics
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
u1.dY = u.dY(:,1);
u1.f = u.f(:,1);
%User Line: u1 = u(:,1);
cada1temp1 = Gator1Data.Data1;
dx.dY = x2.dY;
dx.f = cada1temp1;
dx.f(:,1) = x2.f;
%User Line: dx(:,1) = x2;
cada2f1 = size(dx.f,1);
cada1td1 = zeros(cada2f1,2);
cada1td1(:,2) = u1.dY;
cada2f1 = dx.dY(:,1);
cada1td1(:,1) = cada2f1;
dx.dY = cada1td1;
dx.f(:,2) = u1.f;
%User Line: dx(:,2) = u1;
%User Line: %------------- END OF CODE --------------
dx.dY_size = [2 3];
dx.dY_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_dynamics_YY
ADiGator_dynamics_YY = load('dynamics_YY.mat');
return
end