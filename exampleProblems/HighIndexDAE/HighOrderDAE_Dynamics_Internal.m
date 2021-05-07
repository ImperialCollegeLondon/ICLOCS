function [dx,geqn] = HighOrderDAE_Dynamics_Internal(x,u,p,t,vdat)
%Aly-Chan problem Dynamics - Internal
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
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
%
%------------- BEGIN CODE --------------

x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);x4 = x(:,4);x5 = u(:,1);
u1 = u(:,2);

dx1=x2;
dx2=-x5.*x1-vdat.a.*x2+u1.*x3;
dx3=x4;
dx4=-vdat.g-2*x5.*x3-u1.*x1-vdat.a.*x4;

g1=x1.^2+x3.^2-vdat.L;

dx=[dx1,dx2,dx3,dx4];

geqn=g1;

%------------- END OF CODE --------------