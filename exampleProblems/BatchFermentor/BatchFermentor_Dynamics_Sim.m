
function dx = BatchFermentor_Dynamics_Sim(x,u,p,t,vdat)
% Fed-batch fermentor Dynamics - Simulation
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
x1 = x(:,1);x2 = x(:,2);x3 = x(:,3);x4 = x(:,4);
u1 = u(:,1);

h1 = 0.11*(x3./(0.006*x1+x3));
h2 = 0.0055*(x3./(0.0001+x3.*(1+10*x3)));

dx(:,1) = (h1.*x1-u1.*(x1./500./x4));
dx(:,2) = (h2.*x1-0.01*x2-u1.*(x2./500./x4));
dx(:,3) = (-h1.*x1/0.47-h2.*x1/1.2-x1.*(0.029*x3./(0.0001+x3))+u1./x4.*(1-x3/500));
dx(:,4) = u1/500;

%------------- END OF CODE --------------