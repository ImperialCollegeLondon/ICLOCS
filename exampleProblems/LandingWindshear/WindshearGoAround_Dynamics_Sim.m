function dx = WindshearGoAround_Dynamics_Sim(x,u,p,t,vdat)
%Aircraft go around in the present of wind-shear problem - Dynamics - simulation
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

auxdata = vdat.auxdata;

pos = x(:,1);
h = x(:,2);
v = x(:,3);
fpa = x(:,4);
alpha = u(:,1);

T=ppval(auxdata.beta_poly,t).*(auxdata.A_0+auxdata.A_1.*v+auxdata.A_2.*v.^2);
D=0.5*(auxdata.B_0+auxdata.B_1.*alpha+auxdata.B_2.*alpha.^2).*auxdata.rho.*auxdata.S.*v.^2;
L=0.5*ppval(auxdata.C_L_poly,alpha).*auxdata.rho.*auxdata.S.*v.^2;

posdot=v.*cos(fpa)+ppval(auxdata.W1_poly,pos);
hdot=v.*sin(fpa)+ppval(auxdata.W2_poly,pos).*h/1000;

W1_dot=ppval(auxdata.W1_poly,pos+posdot)-ppval(auxdata.W1_poly,pos);
W2_dot=ppval(auxdata.W2_poly,pos+posdot).*(h+hdot)/1000-ppval(auxdata.W2_poly,pos).*h/1000;
vdot=T./auxdata.m.*cos(alpha+auxdata.delta)-D./auxdata.m-auxdata.g*sin(fpa)-(W1_dot.*cos(fpa)+W2_dot.*sin(fpa));
fpadot=T./auxdata.m./v.*sin(alpha+auxdata.delta)+L./auxdata.m./v-auxdata.g./v.*cos(fpa)+(W1_dot.*sin(fpa)-W2_dot.*cos(fpa))./v;

dx = [posdot, hdot, vdot, fpadot];


%------------- END OF CODE ----------------