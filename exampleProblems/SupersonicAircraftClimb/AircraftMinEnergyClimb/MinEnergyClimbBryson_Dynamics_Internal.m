function [dx] = MinEnergyClimbBryson_Dynamics_Internal(x,u,p,t,vdat)
%Supersonic Aircraft Minimum Fuel Climb Problem - Dynamics - Internal
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

Atomsrho = vdat.Atomsrho;
Atomssos = vdat.Atomssos;
TLT = vdat.T;
mu = vdat.mu;
S = vdat.S;
g0 = vdat.g0;
Isp = vdat.Isp;
Re = vdat.Re;

h = x(:,1);
v = x(:,2);
fpa = x(:,3);
mass = x(:,4);
alpha = u(:,1);

r = h+Re;
rho = ppval(Atomsrho,h);
sos = ppval(Atomssos,h);
Mach = v./sos;

CD0 = ppval(vdat.CDdat,Mach);
Clalpha = ppval(vdat.Clalphadat,Mach);
eta = ppval(vdat.etadat,Mach);

Thrust = interp2(vdat.aa,vdat.mm,TLT,h,Mach,'spline');
CD = CD0 + eta.*Clalpha.*alpha.^2;
CL = Clalpha.*alpha;
q = 0.5.*rho.*v.*v;
D = q.*S.*CD;
L = q.*S.*CL;

hdot = v.*sin(fpa);
vdot = (Thrust.*cos(alpha)-D)./mass - mu.*sin(fpa)./r.^2;
fpadot = (Thrust.*sin(alpha)+L)./(mass.*v)+cos(fpa).*(v./r-mu./(v.*r.^2));
mdot = -Thrust./(g0.*Isp);

dx = [hdot, vdot, fpadot, mdot];
