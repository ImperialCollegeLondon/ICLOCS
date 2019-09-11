function [dx,g_neq] = CarParking_Dynamics_Internal(x,u,p,t,vdat)
% Dynamics for Internal Model
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

auxdata = vdat.auxdata;

posx = x(:,1);
posy = x(:,2);
v = x(:,3);
theta = x(:,4);
phi = x(:,5);
a=u(:,1);
u2=u(:,2);

posx_dot=v.*cos(theta);
posy_dot=v.*sin(theta);
v_dot=a;
theta_dot=v.*tan(phi)./auxdata.l_axes;
phi_dot=u2;

pw=vdat.penalty.values(vdat.penalty.i);

curvature_dot=u2./auxdata.l_axes./(cos(phi)).^2;

A_x=posx+(auxdata.l_axes+auxdata.l_front).*cos(theta)-auxdata.b_width.*sin(theta);
B_x=posx+(auxdata.l_axes+auxdata.l_front).*cos(theta)+auxdata.b_width.*sin(theta);
C_x=posx-auxdata.l_rear.*cos(theta)+auxdata.b_width.*sin(theta);
D_x=posx-auxdata.l_rear.*cos(theta)-auxdata.b_width.*sin(theta);

f_p_A_x=(-tanh(pw*(A_x-5/pw))+tanh(pw*(A_x-auxdata.SL+5/pw)))/2.*auxdata.SW;
f_p_B_x=(-tanh(pw*(B_x-0.1))+tanh(pw*(B_x-auxdata.SL+5/pw)))/2.*auxdata.SW;
f_p_C_x=(-tanh(pw*(C_x-0.1))+tanh(pw*(C_x-auxdata.SL+5/pw)))/2.*auxdata.SW;
f_p_D_x=(-tanh(pw*(D_x-0.1))+tanh(pw*(D_x-auxdata.SL+5/pw)))/2.*auxdata.SW;

A_y=posy+(auxdata.l_axes+auxdata.l_front).*sin(theta)+auxdata.b_width.*cos(theta);
B_y=posy+(auxdata.l_axes+auxdata.l_front).*sin(theta)-auxdata.b_width.*cos(theta);
C_y=posy-auxdata.l_rear.*sin(theta)-auxdata.b_width.*cos(theta);
D_y=posy-auxdata.l_rear.*sin(theta)+auxdata.b_width.*cos(theta);

Ox_p=-posx.*cos(theta)-posy.*sin(theta)-(auxdata.l_axes+auxdata.l_front-auxdata.l_rear)/2;
Oy_p=posx.*sin(theta)-posy.*cos(theta);

Ex_p=-posx.*cos(theta)-posy.*sin(theta)-(auxdata.l_axes+auxdata.l_front-auxdata.l_rear)/2+auxdata.SL*cos(theta);
Ey_p=posx.*sin(theta)-posy.*cos(theta)-auxdata.SL*sin(theta);

Ax_p=(auxdata.l_axes+auxdata.l_front+auxdata.l_rear)/2*ones(size(Ox_p));
Ay_p=auxdata.b_width*ones(size(Oy_p));
Bx_p=(auxdata.l_axes+auxdata.l_front+auxdata.l_rear)/2*ones(size(Ox_p));
By_p=-auxdata.b_width*ones(size(Oy_p));
Cx_p=-(auxdata.l_axes+auxdata.l_front+auxdata.l_rear)/2*ones(size(Ox_p));
Cy_p=-auxdata.b_width*ones(size(Oy_p));
Dx_p=-(auxdata.l_axes+auxdata.l_front+auxdata.l_rear)/2*ones(size(Ox_p));
Dy_p=auxdata.b_width*ones(size(Oy_p));

AreaO1=polyarea([Ox_p,Ax_p,Bx_p],[Oy_p,Ay_p,By_p],2);
AreaO2=polyarea([Ox_p,Cx_p,Bx_p],[Oy_p,Cy_p,By_p],2);
AreaO3=polyarea([Ox_p,Ax_p,Dx_p],[Oy_p,Ay_p,Dy_p],2);
AreaO4=polyarea([Ox_p,Dx_p,Cx_p],[Oy_p,Dy_p,Cy_p],2);
AreaO=AreaO1+AreaO2+AreaO3+AreaO4;

AreaE1=polyarea([Ex_p,Ax_p,Bx_p],[Ey_p,Ay_p,By_p],2);
AreaE2=polyarea([Ex_p,Cx_p,Bx_p],[Ey_p,Cy_p,By_p],2);
AreaE3=polyarea([Ex_p,Ax_p,Dx_p],[Ey_p,Ay_p,Dy_p],2);
AreaE4=polyarea([Ex_p,Dx_p,Cx_p],[Ey_p,Dy_p,Cy_p],2);
AreaE=AreaE1+AreaE2+AreaE3+AreaE4;

AreaRef=(auxdata.l_axes+auxdata.l_front+auxdata.l_rear)*2*auxdata.b_width;

 
%% return dynamics and inequality constraints
dx = [posx_dot, posy_dot, v_dot, theta_dot, phi_dot];
g_neq=[curvature_dot, auxdata.CL-A_y, auxdata.CL-B_y, auxdata.CL-C_y, auxdata.CL-D_y, A_y-f_p_A_x, B_y-f_p_B_x, C_y-f_p_C_x, D_y-f_p_D_x, AreaO-AreaRef, AreaE-AreaRef ];

