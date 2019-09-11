function dx = OrbitRaising_Dynamics_Sim(x,u,p,t,data)
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

r=x(:,1);q=x(:,2);v=x(:,3);phi=u(:,1);

T1=data.T1;
md=data.md;
dx=[q,...
    (v.*v)./r-r.^(-2)+(T1*sin(phi)./(1-md*t)),...
    (-q.*v)./r+(T1*cos(phi)./(1-md*t))];
