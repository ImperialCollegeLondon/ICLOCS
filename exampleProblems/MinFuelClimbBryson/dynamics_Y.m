% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function dx = dynamics_Y(x,u,p,t,data)
global ADiGator_dynamics_Y
if isempty(ADiGator_dynamics_Y); ADiGator_LoadData(); end
Gator1Data = ADiGator_dynamics_Y.dynamics_Y.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Aircraft Dynamics
%User Line: %
%User Line: % Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
%User Line: % The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
%User Line: % This code is published under the BSD License.
%User Line: % Department of Aeronautics and Department of Electrical and Electronic Engineering,
%User Line: % Imperial College London London  England, UK
%User Line: % ICLOCS (Imperial College London Optimal Control) Version 2.0
%User Line: % 1 May 2018
%User Line: % iclocs@imperial.ac.uk
Atomsrho = data.Atomsrho;
%User Line: Atomsrho = data.Atomsrho;
Atomssos = data.Atomssos;
%User Line: Atomssos = data.Atomssos;
TLT = data.T;
%User Line: TLT = data.T;
mu = data.mu;
%User Line: mu = data.mu;
S = data.S;
%User Line: S = data.S;
g0 = data.g0;
%User Line: g0 = data.g0;
Isp = data.Isp;
%User Line: Isp = data.Isp;
Re = data.Re;
%User Line: Re = data.Re;
h.dY = x.dY(:,1);
h.f = x.f(:,1);
%User Line: h = x(:,1);
v.dY = x.dY(:,2);
v.f = x.f(:,2);
%User Line: v = x(:,2);
fpa.dY = x.dY(:,3);
fpa.f = x.f(:,3);
%User Line: fpa = x(:,3);
mass.dY = x.dY(:,4);
mass.f = x.f(:,4);
%User Line: mass = x(:,4);
alpha.dY = u.dY(:,1);
alpha.f = u.f(:,1);
%User Line: alpha = u(:,1);
r.dY = h.dY;
r.f = h.f + Re;
%User Line: r = h+Re;
cada1tf1 = ppval(Gator1Data.Data1,h.f);
rho.dY = cada1tf1(:,1).*h.dY;
rho.f = ppval(Atomsrho,h.f);
%User Line: rho = ppval(Atomsrho,h);
cada1tf1 = ppval(Gator1Data.Data2,h.f);
sos.dY = cada1tf1(:,1).*h.dY;
sos.f = ppval(Atomssos,h.f);
%User Line: sos = ppval(Atomssos,h);
cada1td1 = zeros(size(v.dY,1),2);
cada1td1(:,2) = v.dY./sos.f;
cada1td1(:,1) = cada1td1(:,1) + -v.f./sos.f.^2.*sos.dY;
Mach.dY = cada1td1;
Mach.f = v.f./sos.f;
%User Line: Mach = v./sos;
ii.f = ge(Mach.f,0.8);
%User Line: ii = Mach>=0.8;
jj.f = lt(Mach.f,0.8);
%User Line: jj = Mach<0.8;
mpoly.dY = Mach.dY(:,Gator1Data.Index1);
mpoly.f = Mach.f(:);
%User Line: mpoly = Mach(ii);
cada1f1 = length(Mach.f);
CD0.f = zeros(cada1f1,1);
%User Line: CD0 = zeros(length(Mach),1);
cada1f1 = length(Mach.f);
Clalpha.f = zeros(cada1f1,1);
%User Line: Clalpha = zeros(length(Mach),1);
cada1f1 = length(Mach.f);
eta.f = zeros(cada1f1,1);
%User Line: eta = zeros(length(Mach),1);
cadaconditional1 = any(ii.f,1);
%User Line: cadaconditional1 = any(ii);
if cadaconditional1
    cada1tf1 = ppval(Gator1Data.Data3,mpoly.f);
    cada1f1dY = cada1tf1(:,Gator1Data.Index2).*mpoly.dY;
    cada1f1 = ppval(data.CDdat,mpoly.f);
    cada1td2 = zeros(size(CD0.f,1),2);
    cada1tind1 = ii.f(:,Gator1Data.Index3);
    cada1td2(cada1tind1) = cada1f1dY(cada1tind1);
    CD0.dY = cada1td2;
    CD0.f(ii.f) = cada1f1(ii.f);
    %User Line: CD0(ii) = ppval(data.CDdat,mpoly);
    cada1tf1 = ppval(Gator1Data.Data4,mpoly.f);
    cada1f1dY = cada1tf1(:,Gator1Data.Index4).*mpoly.dY;
    cada1f1 = ppval(data.Clalphadat,mpoly.f);
    cada1td2 = zeros(size(Clalpha.f,1),2);
    cada1tind1 = ii.f(:,Gator1Data.Index5);
    cada1td2(cada1tind1) = cada1f1dY(cada1tind1);
    Clalpha.dY = cada1td2;
    Clalpha.f(ii.f) = cada1f1(ii.f);
    %User Line: Clalpha(ii) = ppval(data.Clalphadat,mpoly);
    cada1tf1 = ppval(Gator1Data.Data5,mpoly.f);
    cada1f1dY = cada1tf1(:,Gator1Data.Index6).*mpoly.dY;
    cada1f1 = ppval(data.etadat,mpoly.f);
    cada1td2 = zeros(size(eta.f,1),2);
    cada1tind1 = ii.f(:,Gator1Data.Index7);
    cada1td2(cada1tind1) = cada1f1dY(cada1tind1);
    eta.dY = cada1td2;
    eta.f(ii.f) = cada1f1(ii.f);
    %User Line: eta(ii) = ppval(data.etadat,mpoly);
else
    CD0.dY = zeros(size(CD0.f,1),2);
    Clalpha.dY = zeros(size(Clalpha.f,1),2);
    eta.dY = zeros(size(eta.f,1),2);
end
cadaconditional1 = any(jj.f,1);
%User Line: cadaconditional1 = any(jj);
if cadaconditional1
    cada1td2 = zeros(size(CD0.f,1),2);
    cada1td2(:,Gator1Data.Index8) = CD0.dY;
    cada1tind1 = jj.f(:,Gator1Data.Index9);
    cada1td2(cada1tind1) = 0;
    CD0.dY = cada1td2;
    CD0.f(jj.f) = 0.013;
    %User Line: CD0(jj) = 0.013;
    cada1td2 = zeros(size(Clalpha.f,1),2);
    cada1td2(:,Gator1Data.Index10) = Clalpha.dY;
    cada1tind1 = jj.f(:,Gator1Data.Index11);
    cada1td2(cada1tind1) = 0;
    Clalpha.dY = cada1td2;
    Clalpha.f(jj.f) = 3.44;
    %User Line: Clalpha(jj) = 3.44;
    cada1td2 = zeros(size(eta.f,1),2);
    cada1td2(:,Gator1Data.Index12) = eta.dY;
    cada1tind1 = jj.f(:,Gator1Data.Index13);
    cada1td2(cada1tind1) = 0;
    eta.dY = cada1td2;
    eta.f(jj.f) = 0.54;
    %User Line: eta(jj) = 0.54;
else
end
cada1tf1 = adigatorEvalInterp2pp(Gator1Data.Data7,h.f,Mach.f);
cada1tf2 = adigatorEvalInterp2pp(Gator1Data.Data8,h.f,Mach.f);
cada1tf3 = cada1tf1(:,1);
cada1td3 = zeros(size(h.dY,1),2);
cada1td3(:,1) = cada1tf3.*h.dY;
cada1tf3 = cada1tf2(:,Gator1Data.Index14);
cada1td3 = cada1td3 + cada1tf3.*Mach.dY;
Thrust.dY = cada1td3;
Thrust.f = adigatorEvalInterp2pp(Gator1Data.Data6,h.f,Mach.f);
%User Line: Thrust = interp2(data.aa,data.mm,TLT,h,Mach,'spline');
cada1tf1 = Clalpha.f(:,Gator1Data.Index15);
cada1td1 = cada1tf1.*eta.dY;
cada1tf1 = eta.f(:,Gator1Data.Index16);
cada1td1 = cada1td1 + cada1tf1.*Clalpha.dY;
cada1f1dY = cada1td1;
cada1f1 = eta.f.*Clalpha.f;
cada1f2dY = 2.*alpha.f.^(2-1).*alpha.dY;
cada1f2 = alpha.f.^2;
cada1tf1 = cada1f2(:,Gator1Data.Index17);
cada1td1 = zeros(size(cada1f1dY,1),3);
cada1td1(:,Gator1Data.Index18) = cada1tf1.*cada1f1dY;
cada1td1(:,3) = cada1td1(:,3) + cada1f1.*cada1f2dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f1.*cada1f2;
cada1td1 = zeros(size(CD0.dY,1),3);
cada1td1(:,Gator1Data.Index19) = CD0.dY;
cada1td1 = cada1td1 + cada1f3dY;
CD.dY = cada1td1;
CD.f = CD0.f + cada1f3;
%User Line: CD = CD0 + eta.*Clalpha.*alpha.^2;
cada1tf1 = alpha.f(:,Gator1Data.Index20);
cada1td1 = zeros(size(Clalpha.dY,1),3);
cada1td1(:,Gator1Data.Index21) = cada1tf1.*Clalpha.dY;
cada1td1(:,3) = cada1td1(:,3) + Clalpha.f.*alpha.dY;
CL.dY = cada1td1;
CL.f = Clalpha.f.*alpha.f;
%User Line: CL = Clalpha.*alpha;
cada1f1dY = 0.5.*rho.dY;
cada1f1 = 0.5*rho.f;
cada1td1 = zeros(size(cada1f1dY,1),2);
cada1td1(:,1) = v.f.*cada1f1dY;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*v.dY;
cada1f2dY = cada1td1;
cada1f2 = cada1f1.*v.f;
cada1tf1 = v.f(:,Gator1Data.Index22);
cada1td1 = cada1tf1.*cada1f2dY;
cada1td1(:,2) = cada1td1(:,2) + cada1f2.*v.dY;
q.dY = cada1td1;
q.f = cada1f2.*v.f;
%User Line: q = 0.5.*rho.*v.*v;
cada1f1dY = S.*q.dY;
cada1f1 = q.f*S;
cada1tf1 = CD.f(:,Gator1Data.Index23);
cada1td1 = zeros(size(cada1f1dY,1),3);
cada1td1(:,Gator1Data.Index24) = cada1tf1.*cada1f1dY;
cada1tf1 = cada1f1(:,Gator1Data.Index25);
cada1td1 = cada1td1 + cada1tf1.*CD.dY;
D.dY = cada1td1;
D.f = cada1f1.*CD.f;
%User Line: D = q.*S.*CD;
cada1f1dY = S.*q.dY;
cada1f1 = q.f*S;
cada1tf1 = CL.f(:,Gator1Data.Index26);
cada1td1 = zeros(size(cada1f1dY,1),3);
cada1td1(:,Gator1Data.Index27) = cada1tf1.*cada1f1dY;
cada1tf1 = cada1f1(:,Gator1Data.Index28);
cada1td1 = cada1td1 + cada1tf1.*CL.dY;
L.dY = cada1td1;
L.f = cada1f1.*CL.f;
%User Line: L = q.*S.*CL;
cada1f1dY = cos(fpa.f).*fpa.dY;
cada1f1 = sin(fpa.f);
cada1td1 = zeros(size(v.dY,1),2);
cada1td1(:,1) = cada1f1.*v.dY;
cada1td1(:,2) = cada1td1(:,2) + v.f.*cada1f1dY;
hdot.dY = cada1td1;
hdot.f = v.f.*cada1f1;
%User Line: hdot = v.*sin(fpa);
cada1f1dY = -sin(alpha.f).*alpha.dY;
cada1f1 = cos(alpha.f);
cada1tf1 = cada1f1(:,Gator1Data.Index29);
cada1td1 = zeros(size(Thrust.dY,1),3);
cada1td1(:,Gator1Data.Index30) = cada1tf1.*Thrust.dY;
cada1td1(:,3) = cada1td1(:,3) + Thrust.f.*cada1f1dY;
cada1f2dY = cada1td1;
cada1f2 = Thrust.f.*cada1f1;
cada1td1 = cada1f2dY;
cada1td1 = cada1td1 + -D.dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f2 - D.f;
cada1tf1 = mass.f(:,Gator1Data.Index31);
cada1td1 = zeros(size(cada1f3dY,1),4);
cada1td1(:,Gator1Data.Index32) = cada1f3dY./cada1tf1;
cada1td1(:,3) = cada1td1(:,3) + -cada1f3./mass.f.^2.*mass.dY;
cada1f4dY = cada1td1;
cada1f4 = cada1f3./mass.f;
cada1f5dY = cos(fpa.f).*fpa.dY;
cada1f5 = sin(fpa.f);
cada1f6dY = mu.*cada1f5dY;
cada1f6 = mu*cada1f5;
cada1f7dY = 2.*r.f.^(2-1).*r.dY;
cada1f7 = r.f.^2;
cada1td1 = zeros(size(cada1f6dY,1),2);
cada1td1(:,2) = cada1f6dY./cada1f7;
cada1td1(:,1) = cada1td1(:,1) + -cada1f6./cada1f7.^2.*cada1f7dY;
cada1f8dY = cada1td1;
cada1f8 = cada1f6./cada1f7;
cada1td1 = zeros(size(cada1f4dY,1),5);
cada1td1(:,Gator1Data.Index33) = cada1f4dY;
cada1td1(:,Gator1Data.Index34) = cada1td1(:,Gator1Data.Index34) + -cada1f8dY;
vdot.dY = cada1td1;
vdot.f = cada1f4 - cada1f8;
%User Line: vdot = (Thrust.*cos(alpha)-D)./mass - mu.*sin(fpa)./r.^2;
cada1f1dY = cos(alpha.f).*alpha.dY;
cada1f1 = sin(alpha.f);
cada1tf1 = cada1f1(:,Gator1Data.Index35);
cada1td1 = zeros(size(Thrust.dY,1),3);
cada1td1(:,Gator1Data.Index36) = cada1tf1.*Thrust.dY;
cada1td1(:,3) = cada1td1(:,3) + Thrust.f.*cada1f1dY;
cada1f2dY = cada1td1;
cada1f2 = Thrust.f.*cada1f1;
cada1td1 = cada1f2dY;
cada1td1 = cada1td1 + L.dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f2 + L.f;
cada1td1 = zeros(size(mass.dY,1),2);
cada1td1(:,2) = v.f.*mass.dY;
cada1td1(:,1) = cada1td1(:,1) + mass.f.*v.dY;
cada1f4dY = cada1td1;
cada1f4 = mass.f.*v.f;
cada1tf1 = cada1f4(:,Gator1Data.Index37);
cada1td1 = zeros(size(cada1f3dY,1),4);
cada1td1(:,Gator1Data.Index38) = cada1f3dY./cada1tf1;
cada1tf1 = cada1f3(:,Gator1Data.Index39);
cada1tf2 = cada1f4(:,Gator1Data.Index40);
cada1td1(:,Gator1Data.Index41) = cada1td1(:,Gator1Data.Index41) + -cada1tf1./cada1tf2.^2.*cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f3./cada1f4;
cada1f6dY = -sin(fpa.f).*fpa.dY;
cada1f6 = cos(fpa.f);
cada1td1 = zeros(size(v.dY,1),2);
cada1td1(:,2) = v.dY./r.f;
cada1td1(:,1) = cada1td1(:,1) + -v.f./r.f.^2.*r.dY;
cada1f7dY = cada1td1;
cada1f7 = v.f./r.f;
cada1f8dY = 2.*r.f.^(2-1).*r.dY;
cada1f8 = r.f.^2;
cada1td1 = zeros(size(v.dY,1),2);
cada1td1(:,2) = cada1f8.*v.dY;
cada1td1(:,1) = cada1td1(:,1) + v.f.*cada1f8dY;
cada1f9dY = cada1td1;
cada1f9 = v.f.*cada1f8;
cada1tf2 = cada1f9(:,Gator1Data.Index42);
cada1f10dY = -mu./cada1tf2.^2.*cada1f9dY;
cada1f10 = mu./cada1f9;
cada1td1 = cada1f7dY;
cada1td1 = cada1td1 + -cada1f10dY;
cada1f11dY = cada1td1;
cada1f11 = cada1f7 - cada1f10;
cada1td1 = zeros(size(cada1f6dY,1),3);
cada1td1(:,3) = cada1f11.*cada1f6dY;
cada1tf1 = cada1f6(:,Gator1Data.Index43);
cada1td1(:,Gator1Data.Index44) = cada1td1(:,Gator1Data.Index44) + cada1tf1.*cada1f11dY;
cada1f12dY = cada1td1;
cada1f12 = cada1f6.*cada1f11;
cada1td1 = zeros(size(cada1f5dY,1),5);
cada1td1(:,Gator1Data.Index45) = cada1f5dY;
cada1td1(:,Gator1Data.Index46) = cada1td1(:,Gator1Data.Index46) + cada1f12dY;
fpadot.dY = cada1td1;
fpadot.f = cada1f5 + cada1f12;
%User Line: fpadot = (Thrust.*sin(alpha)+L)./(mass.*v)+cos(fpa).*(v./r-mu./(v.*r.^2));
cada1f1dY = -Thrust.dY;
cada1f1 = uminus(Thrust.f);
cada1f2 = g0*Isp;
mdot.dY = cada1f1dY./cada1f2;
mdot.f = cada1f1/cada1f2;
%User Line: mdot = -Thrust./(g0.*Isp);
cada1td1 = zeros(size(hdot.f,1),14);
cada1td1(:,Gator1Data.Index47) = hdot.dY;
cada1td1(:,Gator1Data.Index48) = vdot.dY;
cada1td1(:,Gator1Data.Index49) = fpadot.dY;
cada1td1(:,Gator1Data.Index50) = mdot.dY;
dx.dY = cada1td1;
dx.f = [hdot.f vdot.f fpadot.f mdot.f];
%User Line: dx = [hdot, vdot, fpadot, mdot];
dx.dY_size = [4,5];
dx.dY_location = Gator1Data.Index51;
end


function ADiGator_LoadData()
global ADiGator_dynamics_Y
ADiGator_dynamics_Y = load('dynamics_Y.mat');
return
end