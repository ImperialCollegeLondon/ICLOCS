function dx = CommericalFlightProfile_Dynamics_Sim(x,u,p,t,data)
% Aircraft Dynamics - Simulation
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

auxdata = data.auxdata;

H = x(:,1);
V_tas = x(:,4);
gamma = x(:,5);
chi = x(:,6);
W = x(:,7);
alpha = u(:,1);
phi = u(:,2);
throttle = u(:,3);

Temp=auxdata.Ts+H*auxdata.dTdH;
pressure=auxdata.ps*(Temp./auxdata.Ts).^(-auxdata.g/auxdata.dTdH/auxdata.R);
rho=auxdata.rhos*(Temp./auxdata.Ts).^(-(auxdata.g/auxdata.dTdH/auxdata.R+1));
[ V_cas ] = CAS2TAS( auxdata.kappa, pressure, rho, auxdata.ps, auxdata.rhos, V_tas );

ap=auxdata.a1p.*V_cas.^2+auxdata.a2p.*V_cas+auxdata.a3p;
bp=auxdata.b1p.*V_cas.^2+auxdata.b2p.*V_cas+auxdata.b3p;
cp=auxdata.c1p.*V_cas.^2+auxdata.c2p.*V_cas+auxdata.c3p;
Pmax=ap.*H.^2+bp.*H+cp;

P=(Pmax-auxdata.Pidle).*throttle+auxdata.Pidle;
J=60*V_tas/auxdata.nprop/auxdata.Dprop;
ita=-0.13289*J.^6+1.2536*J.^5-4.8906*J.^4+10.146*J.^3-11.918*J.^2+7.6740*J-1.3452;
Thrust=2*745.6*ita.*P./V_tas;

% Calculate aerodynamic forces
cl=auxdata.clalpha.*(alpha-auxdata.alpha0)*180/pi;
L=0.5.*cl.*rho.*V_tas.^2.*auxdata.S;
cd=auxdata.cd0+auxdata.k_cd*cl.^2;
Drag=0.5.*cd.*rho.*V_tas.^2.*auxdata.S;

%% equations of motions
TAS_dot=(Thrust-Drag-W.*sin(gamma)).*auxdata.g./W;
gamma_dot=(L.*cos(phi)+Thrust*sin(auxdata.alphat)-W.*cos(gamma))*auxdata.g./W./V_tas;
H_dot=V_tas.*sin(gamma);
x_dot=V_tas.*cos(gamma).*cos(chi);
y_dot=V_tas.*cos(gamma).*sin(chi);
chi_dot=L.*sin(phi)./cos(gamma)*auxdata.g./W./V_tas;
W_dot= calcWeight(V_cas,H,throttle,auxdata.FFModel);

dx = [H_dot, x_dot, y_dot, TAS_dot, gamma_dot, chi_dot, W_dot];