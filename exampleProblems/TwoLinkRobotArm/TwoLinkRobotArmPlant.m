
function dx = TwoLinkRobotArmPlant(x,u,p,t,vdat)
% Two-link Robot Arm Dynamics
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk
%
%------------- BEGIN CODE --------------
x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);
u1 = u(:,1);u2 = u(:,2);

dx(:,1) = ( sin(x3).*(9/4*cos(x3).*x1.^2)+2*x2.^2 + 4/3*(u1-u2) ...
          - 3/2*cos(x3).*u2 )./ (31/36 + 9/4*sin(x3).^2);

dx(:,2) = -( sin(x3).*(9/4*cos(x3).*x2.^2)+7/2*x1.^2 - 7/3*u2 ...
          + 3/2*cos(x3).*(u1-u2) )./ (31/36 + 9/4*sin(x3).^2);

dx(:,3) = x2-x1;

dx(:,4) = x1;

%------------- END OF CODE --------------