function dx = CarParking_Dynamics_Sim(x,u,p,t,data)
% Dynamics for Simulation
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

auxdata = data.auxdata;

v = x(:,3);
theta = x(:,4);
phi = u(:,2);
a=u(:,1);

posx_dot=v.*cos(theta);
posy_dot=v.*sin(theta);
v_dot=a;
theta_dot=v.*tan(phi)./auxdata.l_axes;


dx = [posx_dot, posy_dot, v_dot, theta_dot];
