function [dx,g_eq] = LTOT_Dynamics_Internal(x,u,p,t,vdat)
%Low Thrust Orbit Transfer Problem - Dynamics - Internal
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
%    vdat - structured variable containing the values of additional vdat used inside
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

% obtain parameters
Isp = vdat.Isp; 
mu = vdat.mu; 
g0 = vdat.g0; 
T = vdat.T; 
Re = vdat.Re; 
J2 = vdat.J2;
J3 = vdat.J3;
J4 = vdat.J4;

tau = p(:,1).*ones(size(x(:,1)));
p = x(:,1);
f = x(:,2);
g = x(:,3);
h = x(:,4);
k = x(:,5);
L = x(:,6);
w = x(:,7);
u_r = u(:,1);
u_t = u(:,2);
u_h = u(:,3);

% define variables as in eqn (6.46-6.41)
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha_sq = h.^2-k.^2;
chi = sqrt(h.^2+k.^2);
s_sq = 1+chi.^2;

% Cartesian state (eqn 6.42-6,43)
r_x = (r./s_sq).*(cos(L)+alpha_sq.*cos(L)+2*h.*k.*sin(L));
r_y = (r./s_sq).*(sin(L)-alpha_sq.*sin(L)+2*h.*k.*cos(L));
r_z = (2*r./s_sq).*(h.*sin(L)-k.*cos(L));
v_x = -(1./s_sq).*sqrt(mu./p).*(sin(L)+alpha_sq.*sin(L)-2*h.*k.*cos(L)+g-2*f.*h.*k+alpha_sq.*g);
v_y = -(1./s_sq).*sqrt(mu./p).*(-cos(L)+alpha_sq.*cos(L)+2*h.*k.*sin(L)-f+2*g.*h.*k+alpha_sq.*f);
v_z = (2./s_sq).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

r_vect = [r_x r_y r_z];
v_vect = [v_x v_y v_z];
r_norm = sqrt(r_x.^2+r_y.^2+r_z.^2);

% cross product for orbital element (eqn 6.45)
rXv = cross(r_vect,v_vect,2);
rXv_norm = sqrt(rXv(:,1).^2+rXv(:,2).^2+rXv(:,3).^2);
rXvXr = cross(rXv,r_vect,2);

i_r1 = r_vect(:,1)./r_norm;
i_r2 = r_vect(:,2)./r_norm;
i_r3 = r_vect(:,3)./r_norm;

i_theta1 = rXvXr(:,1)./(rXv_norm.*r_norm);
i_theta2 = rXvXr(:,2)./(rXv_norm.*r_norm);
i_theta3 = rXvXr(:,3)./(rXv_norm.*r_norm);

i_h1 = rXv(:,1)./rXv_norm;
i_h2 = rXv(:,2)./rXv_norm;
i_h3 = rXv(:,3)./rXv_norm;

i_r = [i_r1 i_r2 i_r3];
i_t = [i_theta1 i_theta2 i_theta3];
i_h = [i_h1 i_h2 i_h3];

% gravitational disturbing acceleration
% eqn 6.47
en_x_ir = i_r(:,3);

en_x_ir_x_ir1 = en_x_ir.*i_r1;
en_x_ir_x_ir2 = en_x_ir.*i_r2;
en_x_ir_x_ir3 = en_x_ir.*i_r3;

en_m_en_x_ir_x_ir1 = 0-en_x_ir_x_ir1;
en_m_en_x_ir_x_ir2 = 0-en_x_ir_x_ir2;
en_m_en_x_ir_x_ir3 = 1-en_x_ir_x_ir3;

en_m_en_x_ir_x_ir_norm = sqrt(en_m_en_x_ir_x_ir1.^2+en_m_en_x_ir_x_ir2.^2+en_m_en_x_ir_x_ir3.^2);
in_1 = en_m_en_x_ir_x_ir1./en_m_en_x_ir_x_ir_norm;
in_2 = en_m_en_x_ir_x_ir2./en_m_en_x_ir_x_ir_norm;
in_3 = en_m_en_x_ir_x_ir3./en_m_en_x_ir_x_ir_norm;

% eqn 6.48-6.49
sinphi = r_z./sqrt(r_x.^2+r_z.^2);
cosphi = sqrt(1-sinphi.^2);

P2 = (3*sinphi.^2-2)./2;
P3 = (5*sinphi.^3-3*sinphi)./2;
P4 = (35*sinphi.^4-30*sinphi.^2+3)./8;
dP2 = 3*sinphi;
dP3 = (15*sinphi-3)./2;
dP4 = (140*sinphi.^3-60*sinphi)./8;

delta_gn = -(mu*cosphi./(r.^2)).*((Re./r).^2.*dP2.*J2+(Re./r).^3.*dP3.*J3+(Re./r).^4.*dP4.*J4);
delta_gr = -(mu./(r.^2)).*((2+1)*(Re./r).^2.*P2.*J2+(3+1)*(Re./r).^3.*P3.*J3+(4+1)*(Re./r).^4.*P4.*J4);

% eqn 6.46
delta_gn_x_in1 = delta_gn.*in_1;
delta_gn_x_in2 = delta_gn.*in_2;
delta_gn_x_in3 = delta_gn.*in_3;
delta_gr_x_ir1 = delta_gr.*i_r1;
delta_gr_x_ir2 = delta_gr.*i_r2;
delta_gr_x_ir3 = delta_gr.*i_r3;
delta_g1 = delta_gn_x_in1 - delta_gr_x_ir1;
delta_g2 = delta_gn_x_in2 - delta_gr_x_ir2;
delta_g3 = delta_gn_x_in3 - delta_gr_x_ir3;

% eqn 6.50
Delta_g1 = i_r(:,1).*delta_g1+i_r(:,2).*delta_g2+i_r(:,3).*delta_g3;
Delta_g2 = i_t(:,1).*delta_g1+i_t(:,2).*delta_g2+i_t(:,3).*delta_g3;
Delta_g3 = i_h(:,1).*delta_g1+i_h(:,2).*delta_g2+i_h(:,3).*delta_g3;


% Thrust acceleration - burn arcs (eqn 6.51)
Delta_T1 = ((g0*T*(1+0.01*tau))./w).*u_r;
Delta_T2 = ((g0*T*(1+0.01*tau))./w).*u_t;
Delta_T3 = ((g0*T*(1+0.01*tau))./w).*u_h;

% Total disburbing acceleration (eqn 6.44)
Delta_1 = Delta_g1+Delta_T1;
Delta_2 = Delta_g2+Delta_T2;
Delta_3 = Delta_g3+Delta_T3;

% equations of motion (eqn 6.31-6.36)
dp = (2*p./q).*sqrt(p./mu).*Delta_2;
df = sqrt(p./mu).*sin(L).*Delta_1+sqrt(p./mu).*(1./q).*((q+1).*cos(L)+f).*Delta_2-sqrt(p./mu).*(g./q).*(h.*sin(L)-k.*cos(L)).*Delta_3;
dg = -sqrt(p./mu).*cos(L).*Delta_1+sqrt(p./mu).*(1./q).*((q+1).*sin(L)+g).*Delta_2+sqrt(p./mu).*(f./q).*(h.*sin(L)-k.*cos(L)).*Delta_3;
dh = sqrt(p./mu).*(s_sq.*cos(L)./(2*q)).*Delta_3;
dk = sqrt(p./mu).*(s_sq.*sin(L)./(2*q)).*Delta_3;
dL = sqrt(p./mu).*(1./q).*(h.*sin(L)-k.*cos(L)).*Delta_3+sqrt(mu.*p).*((q./p).^2);
dw = -(T*(1+0.01*tau)/Isp);

% Return variables
dx=[dp,df,dg,dh,dk,dL,dw];
g_eq=[u_r.^2+u_t.^2+u_h.^2-1];