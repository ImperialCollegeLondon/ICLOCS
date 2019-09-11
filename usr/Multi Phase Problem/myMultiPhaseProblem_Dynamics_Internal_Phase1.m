
function [dx,g_eq,g_neq] = myMultiPhaseProblem_Dynamics_Internal_Phase1(x,u,p,t,vdat)
% Template for specifying the dynamics for internal model 
%
% Syntax:  
%          [dx] = myProblem_Dynamics_Internal(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = myProblem_Dynamics_Internal(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
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
%    g_eq - constraint function for equality constraints
%    g_neq - constraint function for inequality constraints
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

%Define Path constraints
g_eq(:,1)=g_eq1(x1,...,u1,...p,t);
g_eq(:,2)=g_eq2(x1,...,u1,...p,t);
...

g_neq(:,1)=g_neq1(x1,...,u1,...p,t);
g_neq(:,2)=g_neq2(x1,...,u1,...p,t);
...

%------------- END OF CODE --------------