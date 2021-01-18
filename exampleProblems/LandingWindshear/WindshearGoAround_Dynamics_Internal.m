function [dx,g_neq] = WindshearGoAround_Dynamics_Internal(x,u,p,t,vdat)
%Aircraft go around in the present of wind-shear problem - Dynamics - Internal
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

auxdata = vdat.auxdata;

pos = x(:,1);
h = x(:,2);
v = x(:,3);
fpa = x(:,4);
alpha = u(:,1);

T=ppval(auxdata.beta_poly,t).*(auxdata.A_0+auxdata.A_1.*v+auxdata.A_2.*v.^2);
D=0.5*(auxdata.B_0+auxdata.B_1.*alpha+auxdata.B_2.*alpha.^2).*auxdata.rho.*auxdata.S.*v.^2;
L=0.5*ppval(auxdata.C_L_poly,alpha).*auxdata.rho.*auxdata.S.*v.^2;

posdot=v.*cos(fpa)+x(:,5);
hdot=v.*sin(fpa)+x(:,6);

W1_dot=ppval(auxdata.W1_dot_poly,pos).*posdot;
W2_dot=zeros(size(W1_dot));
idx=pos>=500 & pos<=4100;
W2_dot(idx)=(204*auxdata.c.*exp(-auxdata.c.*(pos(idx) - 2300).^4).*(pos(idx) - 2300).^3).*posdot(idx).*h(idx)/1000+(-51.*exp(-auxdata.c.*(pos(idx)-2300).^4))/1000.*hdot(idx);
W2_dot(~idx)=ppval(auxdata.W2_dot_poly,pos(~idx)).*posdot(~idx).*h(~idx)/1000+ppval(auxdata.W2_poly,pos(~idx))/1000.*hdot(~idx);

vdot=T./auxdata.m.*cos(alpha+auxdata.delta)-D./auxdata.m-auxdata.g*sin(fpa)-(W1_dot.*cos(fpa)+W2_dot.*sin(fpa));
fpadot=T./auxdata.m./v.*sin(alpha+auxdata.delta)+L./auxdata.m./v-auxdata.g./v.*cos(fpa)+(W1_dot.*sin(fpa)-W2_dot.*cos(fpa))./v;

dx = [posdot, hdot, vdot, fpadot, W1_dot W2_dot];
g_neq=h-p(1);

%------------- END OF CODE --------------