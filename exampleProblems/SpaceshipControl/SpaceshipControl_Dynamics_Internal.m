function [dx] = SpaceshipControl_Dynamics_Internal(x,u,p,t,data)
% Optimal Control for a Spaceship - Dynamics - Internal
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

px = x(:,1);py = x(:,2);vx = x(:,3);vy = x(:,4);
ux = u(:,1);uy = u(:,2);

mu=data.mu;
x1=-mu;
x2=1-mu;
r1=sqrt((px-x1).^2+py.^2);
r2=sqrt((px-x2).^2+py.^2);

dx(:,1) = vx;

dx(:,2) = vy;

dx(:,3) = 2.*vy+px-(1-mu).*(px-x1)./(r1.^3)-mu.*(px-x2)./(r2.^3)+ux;

dx(:,4) = -2.*vx+py-(1-mu).*py./(r1.^3)-mu.*py./(r2.^3)+uy;

dx(:,5) = ux.^2+uy.^2;
