
function [dx] = CartPoleSwingUp_Dynamics_Internal_Phase2(x,u,p,t,vdat)
% Cart Pole Swing-up Problem Dynamics - Internal
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
x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);x4 = x(:,4);
u1 = u(:,1);


dx(:,1) = (vdat.L.*vdat.m2.*sin(x4).*x2.^2 + u1 + vdat.m2.*vdat.g.*cos(x4).*sin(x4))./...
            (vdat.m1+vdat.m2.*(1-cos(x4).^2));

dx(:,2) = -( vdat.L.*vdat.m2.*cos(x4).*sin(x4).*x2.^2 + u1.*cos(x4) + (vdat.m1+vdat.m2).*vdat.g.*sin(x4))./...
            (vdat.L.*vdat.m1+vdat.L.*vdat.m2.*(1-cos(x4).^2));

dx(:,3) = x1;

dx(:,4) = x2;


%------------- END OF CODE --------------