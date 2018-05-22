
function dx = BangBangPlant(x,u,p,t,vdat)
%Double Integrator Dynamics
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

x1 = x(:,1);x2 = x(:,2);u1 = u(:,1);

dx(:,1) = x2;

dx(:,2) = u1;

%------------- END OF CODE --------------