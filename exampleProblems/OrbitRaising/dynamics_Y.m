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
%User Line: % f - Returns the ODE right hand side where x'= f(x,u,p,t)
%User Line: % The function must be vectorized and
%User Line: % xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
%User Line: % state corresponds to one column of dx.
%User Line: %
%User Line: %
%User Line: % Syntax:  dx = f(x,u,p,t,data)
%User Line: %
%User Line: % Inputs:
%User Line: %    x  - state vector
%User Line: %    u  - input
%User Line: %    p  - parameter
%User Line: %    t  - time
%User Line: %    data-structured variable containing the values of additional data used inside
%User Line: %          the function
%User Line: %
%User Line: % Output:
%User Line: %    dx - time derivative of x
%User Line: %
%User Line: %  Remark: If the i-th ODE right hand side does not depend on variables it is necessary to multiply
%User Line: %          the assigned value by a vector of ones with the same length  of t  in order
%User Line: %          to have  a vector with the right dimesion  when called for the optimization.
%User Line: %          Example: dx(:,i)= 0*ones(size(t,1));
%User Line: %
%User Line: %------------- BEGIN CODE --------------
r.dY = x.dY(:,1);
r.f = x.f(:,1);
%User Line: r=x(:,1);
theta.dY = x.dY(:,2);
theta.f = x.f(:,2);
%User Line: theta=x(:,2);
v_r.dY = x.dY(:,3);
v_r.f = x.f(:,3);
%User Line: v_r=x(:,3);
v_theta.dY = x.dY(:,4);
v_theta.f = x.f(:,4);
%User Line: v_theta=x(:,4);
u_1.dY = u.dY(:,1);
u_1.f = u.f(:,1);
%User Line: u_1=u(:,1);
u_2.dY = u.dY(:,2);
u_2.f = u.f(:,2);
%User Line: u_2=u(:,2);
T1 = data.T1;
%User Line: T1=data.T1;
md = data.md;
%User Line: md=data.md;
cada1td1 = zeros(size(v_theta.dY,1),2);
cada1td1(:,2) = v_theta.dY./r.f;
cada1td1(:,1) = cada1td1(:,1) + -v_theta.f./r.f.^2.*r.dY;
cada1f1dY = cada1td1;
cada1f1 = v_theta.f./r.f;
cada1td1 = v_theta.f.*v_theta.dY;
cada1td1 = cada1td1 + v_theta.f.*v_theta.dY;
cada1f2dY = cada1td1;
cada1f2 = v_theta.f.*v_theta.f;
cada1td1 = zeros(size(cada1f2dY,1),2);
cada1td1(:,2) = cada1f2dY./r.f;
cada1td1(:,1) = cada1td1(:,1) + -cada1f2./r.f.^2.*r.dY;
cada1f3dY = cada1td1;
cada1f3 = cada1f2./r.f;
cada1f4dY = (-2).*r.f.^((-2)-1).*r.dY;
cada1f4dY((r.f == 0 & r.dY == 0) | (-2) == 0) = 0;
cada1f4 = r.f.^(-2);
cada1td1 = cada1f3dY;
cada1td1(:,1) = cada1td1(:,1) + -cada1f4dY;
cada1f5dY = cada1td1;
cada1f5 = cada1f3 - cada1f4;
cada1f6dT = md.*t.dT;
cada1f6 = md*t.f;
cada1f7dT = -cada1f6dT;
cada1f7 = 1 - cada1f6;
cada1f8dT = -T1./cada1f7.^2.*cada1f7dT;
cada1f8 = T1./cada1f7;
cada1f9dY = cada1f8.*u_1.dY;
cada1f9dT = u_1.f.*cada1f8dT;
cada1f9 = cada1f8.*u_1.f;
cada1td1 = zeros(size(cada1f5dY,1),3);
cada1td1(:,Gator1Data.Index1) = cada1f5dY;
cada1td1(:,3) = cada1td1(:,3) + cada1f9dY;
cada1f10dY = cada1td1;
cada1f10dT = cada1f9dT;
cada1f10 = cada1f5 + cada1f9;
cada1f11dY = -v_r.dY;
cada1f11 = uminus(v_r.f);
cada1td1 = zeros(size(cada1f11dY,1),2);
cada1td1(:,1) = v_theta.f.*cada1f11dY;
cada1td1(:,2) = cada1td1(:,2) + cada1f11.*v_theta.dY;
cada1f12dY = cada1td1;
cada1f12 = cada1f11.*v_theta.f;
cada1tf1 = r.f(:,Gator1Data.Index2);
cada1td1 = zeros(size(cada1f12dY,1),3);
cada1td1(:,Gator1Data.Index3) = cada1f12dY./cada1tf1;
cada1td1(:,1) = cada1td1(:,1) + -cada1f12./r.f.^2.*r.dY;
cada1f13dY = cada1td1;
cada1f13 = cada1f12./r.f;
cada1f14dT = md.*t.dT;
cada1f14 = md*t.f;
cada1f15dT = -cada1f14dT;
cada1f15 = 1 - cada1f14;
cada1f16dT = -T1./cada1f15.^2.*cada1f15dT;
cada1f16 = T1./cada1f15;
cada1f17dY = cada1f16.*u_2.dY;
cada1f17dT = u_2.f.*cada1f16dT;
cada1f17 = cada1f16.*u_2.f;
cada1td1 = zeros(size(cada1f13dY,1),4);
cada1td1(:,Gator1Data.Index4) = cada1f13dY;
cada1td1(:,4) = cada1td1(:,4) + cada1f17dY;
cada1f18dY = cada1td1;
cada1f18dT = cada1f17dT;
cada1f18 = cada1f13 + cada1f17;
cada1td1 = zeros(size(v_r.f,1),10);
cada1td1(:,4) = v_r.dY;
cada1td1(:,Gator1Data.Index5) = cada1f1dY;
cada1td1(:,Gator1Data.Index6) = cada1f10dY;
cada1td1(:,Gator1Data.Index7) = cada1f18dY;
dx.dY = cada1td1;
cada1td1 = zeros(size(v_r.f,1),2);
cada1td1(:,1) = cada1f10dT;
cada1td1(:,2) = cada1f18dT;
dx.dT = cada1td1;
dx.f = [v_r.f cada1f1 cada1f10 cada1f18];
%User Line: dx=[v_r,v_theta./r,(v_theta.*v_theta)./r-r.^(-2)+(T1./(1-md*t)).*u_1,    (-v_r.*v_theta)./r+(T1./(1-md*t)).*u_2];
dx.dY_size = [4,6];
dx.dY_location = Gator1Data.Index8;
dx.dT_size = 4;
dx.dT_location = Gator1Data.Index9;
end


function ADiGator_LoadData()
global ADiGator_dynamics_Y
ADiGator_dynamics_Y = load('dynamics_Y.mat');
return
end