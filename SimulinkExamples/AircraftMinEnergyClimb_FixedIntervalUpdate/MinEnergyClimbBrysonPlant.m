function dx = MinEnergyClimbBrysonPlant(x,u,p,t,data)
% Aircraft Dynamics
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

Atomsrho = data.Atomsrho;
Atomssos = data.Atomssos;
TLT = data.T;
mu = data.mu;
S = data.S;
g0 = data.g0;
Isp = data.Isp;
Re = data.Re;

h = x(:,1);
v = x(:,2);
fpa = x(:,3);
mass = x(:,4);
alpha = u(:,1);

r = h+Re;
rho = ppval(Atomsrho,h);
sos = ppval(Atomssos,h);
Mach = v./sos;

ii = Mach>=0.8;
jj = Mach<0.8;
mpoly = Mach(ii);
CD0 = zeros(length(Mach),1);
Clalpha = zeros(length(Mach),1);
eta = zeros(length(Mach),1);
if any(ii)
CD0(ii) = ppval(data.CDdat,mpoly);
Clalpha(ii) = ppval(data.Clalphadat,mpoly);
eta(ii) = ppval(data.etadat,mpoly);
end
if any(jj)
CD0(jj) = 0.013;
Clalpha(jj) = 3.44;
eta(jj) = 0.54;
end

Thrust = interp2(data.aa,data.mm,TLT,h,Mach,'spline');
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
