
function [dx] = myMultiPhaseProblem_Dynamics_Sim_Phase2(x,u,p,t,vdat)
% Template for specifying the dynamics for simulation 
%
% Syntax:  
%          [dx] = myProblem_Dynamics_Sim(x,u,p,t,vdat)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

%Stored data
auxdata = data.auxdata;

%Define states
x1 = x(:,1);
%...
xn=x(:,n),

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);


%Define ODE right-hand side
dx(:,1) = f1(x1,..xn,u1,..um,p,t);
%...
dx(:,n) = fn(x1,..xn,u1,..um,p,t);


%------------- END OF CODE --------------