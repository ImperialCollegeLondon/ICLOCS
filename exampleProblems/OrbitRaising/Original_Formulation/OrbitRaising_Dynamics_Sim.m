function dx = OrbitRaising_Dynamics_Sim(x,u,p,t,vdat)
% Orbit Raising Problem - Dynamics - simulation
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)
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

r=x(:,1);theta=x(:,2);v_r=x(:,3);v_theta=x(:,4);
u_1=u(:,1);u_2=u(:,2);

T1=vdat.T1;
md=vdat.md;
dx=[v_r,v_theta./r,(v_theta.*v_theta)./r-r.^(-2)+(T1./(1-md*t)).*u_1,...
    (-v_r.*v_theta)./r+(T1./(1-md*t)).*u_2];
