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

function dx = dynamics_YY(x,u,p,t,data)
global ADiGator_dynamics_YY
if isempty(ADiGator_dynamics_YY); ADiGator_LoadData(); end
Gator1Data = ADiGator_dynamics_YY.dynamics_YY.Gator1Data;
Gator2Data = ADiGator_dynamics_YY.dynamics_YY.Gator2Data;
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
cada2f1 = size(v_theta.dY,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dY = -v_theta.dY./r.f.^2.*r.dY;
cada2f1 = v_theta.dY./r.f;
cada1td1dY = cada2f1dY;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dY = -v_theta.dY;
cada2f2 = uminus(v_theta.f);
cada2f3dY = 2.*r.f.^(2-1).*r.dY;
cada2f3 = r.f.^2;
cada2td1 = zeros(size(cada2f2dY,1),2);
cada2td1(:,2) = cada2f2dY./cada2f3;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dY;
cada2f4dY = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = r.dY(:,Gator2Data.Index1);
cada2f5dY = cada2tf1.*cada2f4dY;
cada2f5 = cada2f4.*r.dY;
cada2f6dY = cada2f5dY;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index2) = cada2f6dY;
cada2td1(:,2) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f1dYdY = cada1td1dY; cada1f1dY = cada1td1;
cada1f1 = v_theta.f./r.f;
cada1td1dY = v_theta.dY.*v_theta.dY;
cada1td1 = v_theta.f.*v_theta.dY;
cada2f1dY = v_theta.dY.*v_theta.dY;
cada2f1 = v_theta.f.*v_theta.dY;
cada2td1 = cada1td1dY;
cada2td1 = cada2td1 + cada2f1dY;
cada1td1dY = cada2td1;
cada1td1 = cada1td1 + cada2f1;
cada1f2dYdY = cada1td1dY; cada1f2dY = cada1td1;
cada1f2 = v_theta.f.*v_theta.f;
cada2f1 = size(cada1f2dY,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f2dYdY,1),2);
cada2td1(:,2) = cada1f2dYdY./r.f;
cada2td1(:,1) = cada2td1(:,1) + -cada1f2dY./r.f.^2.*r.dY;
cada2f1dY = cada2td1;
cada2f1 = cada1f2dY./r.f;
cada1td1dY = cada2f1dY;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dY = -cada1f2dY;
cada2f2 = uminus(cada1f2);
cada2f3dY = 2.*r.f.^(2-1).*r.dY;
cada2f3 = r.f.^2;
cada2td1 = zeros(size(cada2f2dY,1),2);
cada2td1(:,2) = cada2f2dY./cada2f3;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dY;
cada2f4dY = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = r.dY(:,Gator2Data.Index3);
cada2f5dY = cada2tf1.*cada2f4dY;
cada2f5 = cada2f4.*r.dY;
cada2f6dY = cada2f5dY;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index4) = cada2f6dY;
cada2td1(:,Gator2Data.Index5) = cada1td1dY(:,Gator2Data.Index6);
cada1td1dY = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f3dYdY = cada1td1dY; cada1f3dY = cada1td1;
cada1f3 = cada1f2./r.f;
cada2f1dY = (-3).*r.f.^((-3)-1).*r.dY;
cada2f1dY((r.f == 0 & r.dY == 0) | (-3) == 0) = 0;
cada2f1 = r.f.^(-3);
cada2f2dY = (-2).*cada2f1dY;
cada2f2 = (-2)*cada2f1;
cada1f4dYdY = r.dY.*cada2f2dY;
cada1f4dY = cada2f2.*r.dY;
cada2f1 = eq(r.f,0);
cada2f2 = eq(r.dY,0);
cada2f3 = and(cada2f1,cada2f2);
cada2f4 = or(cada2f3,0);
cada2td2 = cada1f4dYdY;
cada2tind1 = cada2f4(:,1);
cada2td2(cada2tind1) = 0;
cada1f4dYdY = cada2td2;
cada1f4dY(cada2f4) = 0;
cada1f4 = r.f.^(-2);
cada1td1dY = cada1f3dYdY; cada1td1 = cada1f3dY;
cada2f1dY = cada1td1dY(:,Gator2Data.Index7);
cada2f1 = cada1td1(:,1);
cada2f2dY = -cada1f4dYdY;
cada2f2 = uminus(cada1f4dY);
cada2td1 = cada2f1dY;
cada2td1(:,1) = cada2td1(:,1) + cada2f2dY;
cada2f3dY = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index8) = cada2f3dY;
cada2td1(:,Gator2Data.Index9) = cada1td1dY(:,Gator2Data.Index10);
cada1td1dY = cada2td1;
cada1td1(:,1) = cada2f3;
cada1f5dYdY = cada1td1dY; cada1f5dY = cada1td1;
cada1f5 = cada1f3 - cada1f4;
cada1f6dT = md*t.dT;
cada1f6 = md*t.f;
cada1f7dT = uminus(cada1f6dT);
cada1f7 = 1 - cada1f6;
cada2f1 = uminus(T1);
cada2f2dT = 2.*cada1f7.^(2-1).*cada1f7dT;
cada2f2 = cada1f7.^2;
cada2f3dT = -cada2f1./cada2f2.^2.*cada2f2dT;
cada2f3 = cada2f1./cada2f2;
cada1f8dTdT = cada1f7dT.*cada2f3dT;
cada1f8dT = cada2f3.*cada1f7dT;
cada1f8 = T1./cada1f7;
cada1f9dYdT = u_1.dY.*cada1f8dT;
cada1f9dY = cada1f8.*u_1.dY;
cada1f9dTdY = cada1f8dT.*u_1.dY;
cada1f9dTdT = u_1.f.*cada1f8dTdT;
cada1f9dT = u_1.f.*cada1f8dT;
cada1f9 = cada1f8.*u_1.f;
cada2f1 = size(cada1f5dY,1);
cada1td1 = zeros(cada2f1,3);
cada1td1dY = cada1f5dYdY;
cada1td1(:,Gator1Data.Index1) = cada1f5dY;
cada2f1 = cada1td1(:,3);
cada2f2dT = cada1f9dYdT;
cada2f2 = cada2f1 + cada1f9dY;
cada1td1dT = cada2f2dT;
cada1td1(:,3) = cada2f2;
cada1f10dYdY = cada1td1dY; cada1f10dYdT = cada1td1dT; cada1f10dY = cada1td1;
cada1f10dTdY = cada1f9dTdY; cada1f10dTdT = cada1f9dTdT; cada1f10dT = cada1f9dT;
cada1f10 = cada1f5 + cada1f9;
cada1f11dY = uminus(v_r.dY);
cada1f11 = uminus(v_r.f);
cada2f1 = size(cada1f11dY,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dY = cada1f11dY.*v_theta.dY;
cada2f1 = v_theta.f.*cada1f11dY;
cada1td1dY = cada2f1dY;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2f2dY = v_theta.dY.*cada1f11dY;
cada2f2 = cada1f11.*v_theta.dY;
cada2f3dY = cada2f2dY;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),2);
cada2td1(:,1) = cada2f3dY;
cada2td1(:,2) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f12dYdY = cada1td1dY; cada1f12dY = cada1td1;
cada1f12 = cada1f11.*v_theta.f;
cada1tf1dY = r.dY(:,Gator2Data.Index11);
cada1tf1 = r.f(:,Gator1Data.Index2);
cada2f1 = size(cada1f12dY,1);
cada1td1 = zeros(cada2f1,3);
cada2tf1 = cada1tf1(:,Gator2Data.Index12);
cada2td1 = zeros(size(cada1f12dYdY,1),4);
cada2td1(:,Gator2Data.Index13) = cada1f12dYdY./cada2tf1;
cada2td1(:,Gator2Data.Index14) = cada2td1(:,Gator2Data.Index14) + -cada1f12dY./cada1tf1.^2.*cada1tf1dY;
cada2f1dY = cada2td1;
cada2f1 = cada1f12dY./cada1tf1;
cada1td1dY = cada2f1dY;
cada1td1(:,Gator1Data.Index3) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dY = -cada1f12dY;
cada2f2 = uminus(cada1f12);
cada2f3dY = 2.*r.f.^(2-1).*r.dY;
cada2f3 = r.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index15);
cada2td1 = zeros(size(cada2f2dY,1),3);
cada2td1(:,Gator2Data.Index16) = cada2f2dY./cada2tf1;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dY;
cada2f4dY = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = r.dY(:,Gator2Data.Index17);
cada2f5dY = cada2tf1.*cada2f4dY;
cada2f5 = cada2f4.*r.dY;
cada2f6dY = cada2f5dY;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),7);
cada2td1(:,Gator2Data.Index18) = cada2f6dY;
cada2td1(:,Gator2Data.Index19) = cada1td1dY(:,Gator2Data.Index20);
cada1td1dY = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f13dYdY = cada1td1dY; cada1f13dY = cada1td1;
cada1f13 = cada1f12./r.f;
cada1f14dT = md*t.dT;
cada1f14 = md*t.f;
cada1f15dT = uminus(cada1f14dT);
cada1f15 = 1 - cada1f14;
cada2f1 = uminus(T1);
cada2f2dT = 2.*cada1f15.^(2-1).*cada1f15dT;
cada2f2 = cada1f15.^2;
cada2f3dT = -cada2f1./cada2f2.^2.*cada2f2dT;
cada2f3 = cada2f1./cada2f2;
cada1f16dTdT = cada1f15dT.*cada2f3dT;
cada1f16dT = cada2f3.*cada1f15dT;
cada1f16 = T1./cada1f15;
cada1f17dYdT = u_2.dY.*cada1f16dT;
cada1f17dY = cada1f16.*u_2.dY;
cada1f17dTdY = cada1f16dT.*u_2.dY;
cada1f17dTdT = u_2.f.*cada1f16dTdT;
cada1f17dT = u_2.f.*cada1f16dT;
cada1f17 = cada1f16.*u_2.f;
cada2f1 = size(cada1f13dY,1);
cada1td1 = zeros(cada2f1,4);
cada1td1dY = cada1f13dYdY;
cada1td1(:,Gator1Data.Index4) = cada1f13dY;
cada2f1 = cada1td1(:,4);
cada2f2dT = cada1f17dYdT;
cada2f2 = cada2f1 + cada1f17dY;
cada1td1dT = cada2f2dT;
cada1td1(:,4) = cada2f2;
cada1f18dYdY = cada1td1dY; cada1f18dYdT = cada1td1dT; cada1f18dY = cada1td1;
cada1f18dTdY = cada1f17dTdY; cada1f18dTdT = cada1f17dTdT; cada1f18dT = cada1f17dT;
cada1f18 = cada1f13 + cada1f17;
cada2f1 = size(v_r.f,1);
cada1td1 = zeros(cada2f1,10);
cada1td1(:,4) = v_r.dY;
cada1td1dY = cada1f1dYdY;
cada1td1(:,Gator1Data.Index5) = cada1f1dY;
cada2td1 = zeros(size(cada1td1,1),7);
cada2td1(:,Gator2Data.Index21) = cada1f10dYdY;
cada2td1(:,Gator2Data.Index22) = cada1td1dY(:,Gator2Data.Index23);
cada1td1dY = cada2td1;
cada1td1dT = cada1f10dYdT;
cada1td1(:,Gator1Data.Index6) = cada1f10dY;
cada2td1 = zeros(size(cada1td1,1),14);
cada2td1(:,Gator2Data.Index24) = cada1f18dYdY;
cada2td1(:,Gator2Data.Index25) = cada1td1dY(:,Gator2Data.Index26);
cada1td1dY = cada2td1;
cada2td1 = zeros(size(cada1td1,1),2);
cada2td1(:,2) = cada1f18dYdT;
cada2td1(:,1) = cada1td1dT(:,1);
cada1td1dT = cada2td1;
cada1td1(:,Gator1Data.Index7) = cada1f18dY;
dx.dYdY = cada1td1dY; dx.dYdT = cada1td1dT; dx.dY = cada1td1;
cada2f1 = size(v_r.f,1);
cada1td1 = zeros(cada2f1,2);
cada1td1dY = cada1f10dTdY;
cada1td1dT = cada1f10dTdT;
cada1td1(:,1) = cada1f10dT;
cada2td1 = zeros(size(cada1td1,1),2);
cada2td1(:,2) = cada1f18dTdY;
cada2td1(:,1) = cada1td1dY(:,1);
cada1td1dY = cada2td1;
cada2td1 = zeros(size(cada1td1,1),2);
cada2td1(:,2) = cada1f18dTdT;
cada2td1(:,1) = cada1td1dT(:,1);
cada1td1dT = cada2td1;
cada1td1(:,2) = cada1f18dT;
dx.dTdY = cada1td1dY; dx.dTdT = cada1td1dT; dx.dT = cada1td1;
dx.f = [v_r.f cada1f1 cada1f10 cada1f18];
%User Line: dx=[v_r,v_theta./r,(v_theta.*v_theta)./r-r.^(-2)+(T1./(1-md*t)).*u_1,    (-v_r.*v_theta)./r+(T1./(1-md*t)).*u_2];
dx.dY_size = [4 6];
dx.dY_location = Gator1Data.Index8;
dx.dT_size = 4;
dx.dT_location = Gator1Data.Index9;
dx.dTdY_size = [dx.dT_size,6];
dx.dTdY_location = [dx.dT_location(Gator2Data.Index27,:), Gator2Data.Index28];
dx.dTdT_size = dx.dT_size;
dx.dTdT_location = dx.dT_location(Gator2Data.Index29,:);
dx.dYdY_size = [dx.dY_size,6];
dx.dYdY_location = [dx.dY_location(Gator2Data.Index30,:), Gator2Data.Index31];
dx.dYdT_size = dx.dY_size;
dx.dYdT_location = dx.dY_location(Gator2Data.Index32,:);
end


function ADiGator_LoadData()
global ADiGator_dynamics_YY
ADiGator_dynamics_YY = load('dynamics_YY.mat');
return
end