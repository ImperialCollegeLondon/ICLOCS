function dx = ShuttleReentryTrajectory_Dynamics_Sim(x,u,p,t,data)
%Space Shuttle Reentry Trajectory Problem - Dynamics - simulation
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
%
%------------- BEGIN CODE --------------

h = x(:,1);
phi = x(:,2);
theta = x(:,3);
v = x(:,4);
gamma = x(:,5);
psi= x(:,6);
alpha= u(:,1);
beta= u(:,2);

rho = data.rho0.*exp(-h./data.hr);
C_L = data.a0+data.a1.*(alpha*180/pi);
C_D = data.b0+data.b1.*(alpha*180/pi)+data.b2.*(alpha*180/pi).^2;
D = 0.5*rho.*v.^2.*data.S.*C_D;
L = 0.5*rho.*v.^2.*data.S.*C_L;
g = data.mu./((h+data.Re).^2);

dh = v.*sin(gamma);
dphi = v.*cos(gamma).*sin(psi)./((h+data.Re).*cos(theta));
dtheta = v.*cos(gamma).*cos(psi)./(h+data.Re);
dv = -D./data.mass-g.*sin(gamma);
dgamma = L./data.mass./v.*cos(beta)+cos(gamma).*(v./(h+data.Re)-g./v);
dpsi = L.*sin(beta)./cos(gamma)./v./data.mass+v.*cos(gamma).*sin(psi).*tan(theta)./(h+data.Re);

dx = [dh, dphi, dtheta, dv, dgamma, dpsi];