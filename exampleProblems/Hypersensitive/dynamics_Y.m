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
%User Line: %Hypersensitive problem Dynamics
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
cada1f1dY = 3.*x.f.^(3-1).*x.dY;
cada1f1 = x.f.^3;
cada1f2dY = -cada1f1dY;
cada1f2 = uminus(cada1f1);
cada1td1 = zeros(size(cada1f2dY,1),2);
cada1td1(:,1) = cada1f2dY;
cada1td1(:,2) = cada1td1(:,2) + u.dY;
dx.dY = cada1td1;
dx.f = cada1f2 + u.f;
%User Line: dx=-x.^3+u;
%User Line: %------------- END OF CODE --------------
dx.dY_size = 2;
dx.dY_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_dynamics_Y
ADiGator_dynamics_Y = load('dynamics_Y.mat');
return
end