function dx = CarParkingPlant(x,u,p,t,data)
% Vehicle Dynamics
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

auxdata = data.auxdata;

v = x(:,3);
theta = x(:,4);
phi = x(:,5);
a=u(:,1);
u2=u(:,2);

posx_dot=v.*cos(theta);
posy_dot=v.*sin(theta);
v_dot=a;
theta_dot=v.*tan(phi)./auxdata.l_axes;
phi_dot=u2;


dx = [posx_dot, posy_dot, v_dot, theta_dot, phi_dot];

