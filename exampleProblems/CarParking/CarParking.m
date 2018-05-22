function [problem,guess] = CarParking
%BangBang - BangBang Control (Double Integrator Minimum Time Repositioning) Problem
%
% The problem was adapted from 
% B. Li, K. Wang, and Z. Shao, "Time-optimal maneuver planning in automatic parallel parking using a simultaneous dynamic optimization approach". IEEE Transactions on Intelligent Transportation Systems, 17(11), pp.3263-3274, 2016.
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% MAT-files required: none
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
% Plant model name
problem.data.plantmodel = 'CarParkingPlant';

% Scenario Parameters
SL=5;
SW=2;
CL=3.5;
l_front=0.8; 
l_axes=2.5;
l_rear=0.7; 
b_width=1.771/2;
phi_max=deg2rad(33);
a_max=0.75;
v_max=2;
u1_max=0.5;
curvature_dot_max=0.6;

% Store data
auxdata.SL=SL;
auxdata.SW=SW;
auxdata.CL=CL;
auxdata.l_front=l_front;
auxdata.l_axes=l_axes;
auxdata.l_rear=l_rear;
auxdata.b_width=b_width;

% Boundary Conditions 
posx0 = SL+l_rear;
posy0 = 1.5;
theta0=0;
v0=0; % Initial velocity (m/s)
a0=0; % Initial accelration (m/s^2)
phi0 = deg2rad(0); % Initial steering angle (rad)

% Limits on Variables
xmin = -10; xmax = 15;
ymin = -SW; ymax = CL;
vmin = -v_max; vmax = v_max;
amin = -a_max; amax = a_max;
thetamin = -inf; thetamax = inf;
phimin = -phi_max; phimax = phi_max;

%%
%Initial Time. t0<tf
problem.time.t0=0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=10;     
problem.time.tf_max=20; 
guess.tf=14;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl = [];
problem.parameters.pu = [];
guess.parameters = [];

% Initial conditions for system.
problem.states.x0=[posx0 posy0 v0 theta0 phi0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[posx0 posy0 v0 theta0 phi0]; 
problem.states.x0u=[posx0 posy0 v0 theta0 phi0]; 

% State bounds. xl=< x <=xu
problem.states.xl=[xmin ymin vmin thetamin phimin]; 
problem.states.xu=[xmax ymax vmax thetamax phimax]; 

% State rate bounds. xrl=< x <=xru
problem.states.xrl=[-inf -inf -inf -inf -inf]; 
problem.states.xru=[inf inf inf inf inf]; 

% State error bounds
problem.states.xErrorTol=[0.1 0.1 0.1 deg2rad(0.5) deg2rad(0.5)];

% State constraint error bounds
problem.states.xConstraintTol=[0.1 0.1 0.1 deg2rad(0.5) deg2rad(0.5)];
problem.states.xrConstraintTol=[0.1 0.1 0.1 deg2rad(0.5) deg2rad(0.5)];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[l_rear ymin 0 theta0-deg2rad(5) -deg2rad(5)]; 
problem.states.xfu=[SL ymax 0 theta0+deg2rad(5) deg2rad(5)];

% Guess the state trajectories with [x0 xf]
guess.time=[0 guess.tf/2 guess.tf];

guess.states(:,1)=[posx0 l_rear SL-l_axes-l_front];
guess.states(:,2)=[posy0 -SW/2 -SW/2];
guess.states(:,3)=[v0 0 0];
guess.states(:,4)=[theta0 deg2rad(1) 0];
guess.states(:,5)=[phi0 0 0];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[amin -curvature_dot_max*l_axes*cos(phimax)^2];
problem.inputs.uu=[amax curvature_dot_max*l_axes*cos(phimax)^2];

problem.inputs.u0l=[amin -curvature_dot_max*l_axes*cos(phimax)^2];
problem.inputs.u0u=[amax curvature_dot_max*l_axes*cos(phimax)^2];

% Input rate bounds
problem.inputs.url=[-u1_max -inf];
problem.inputs.uru=[u1_max inf];

% Input constraint error bounds
problem.inputs.uConstraintTol=[0.1 deg2rad(0.5)];
problem.inputs.urConstraintTol=[0.1 deg2rad(0.5)];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[a0 0 0];
guess.inputs(:,2)=[0 0 0];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[-curvature_dot_max 0 0 0 0 0 0 0 0 0 0];
problem.constraints.gu=[curvature_dot_max inf inf inf inf inf inf inf inf inf inf];
problem.constraints.gTol=[deg2rad(0.01) 1e-04 1e-04 1e-04 1e-04 1e-01 1e-01 1e-01 1e-01 1e-04 1e-04];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[-inf, -inf, -inf, -inf];
problem.constraints.bu=[0 0 0 0];

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;

% For algebraic variable rate constraint
problem.data.xrl=problem.states.xrl;
problem.data.xru=problem.states.xru;
problem.data.xrConstraintTol=problem.states.xrConstraintTol;
problem.data.url=problem.inputs.url;
problem.data.uru=problem.inputs.uru;
problem.data.urConstraintTol=problem.inputs.urConstraintTol;

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc_unscaled,@b_unscaled};
problem.constraintErrorTol=[problem.constraints.gTol,problem.constraints.gTol,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];
%------------- END OF CODE --------------

function stageCost=L_unscaled(x,u,p,t,data)


% L_unscaled - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;
%
%------------- BEGIN CODE --------------


stageCost = 0*t;

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 

% E_unscaled - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

boundaryCost=tf;

%------------- END OF CODE --------------


function dx = f_unscaled(x,u,p,t,data)

% f_unscaled - Returns the ODE right hand side where x'= f(x,u,p,t)
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = f(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    data-structured variable containing the values of additional data used inside
%          the function 
%
% Output:
%    dx - time derivative of x
%
%  Remark: If the i-th ODE right hand side does not depend on variables it is necessary to multiply
%          the assigned value by a vector of ones with the same length  of t  in order 
%          to have  a vector with the right dimesion  when called for the optimization. 
%          Example: dx(:,i)= 0*ones(size(t,1)); 
%
%------------- BEGIN CODE --------------

auxdata = data.auxdata;

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


dx = [posx_dot, posy_dot, v_dot, theta_dot, phi_dot];

%------------- END OF CODE --------------


function c=g_unscaled(x,u,p,t,data)

% g_unscaled - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  c=g(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%   data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    c - constraint function
%
%
%------------- BEGIN CODE --------------
auxdata = data.auxdata;

posx = x(:,1);
posy = x(:,2);
theta = x(:,4);
phi = x(:,5);
u2=u(:,2);

curvature_dot=u2./auxdata.l_axes./(cos(phi)).^2;

A_x=posx+(auxdata.l_axes+auxdata.l_front).*cos(theta)-auxdata.b_width.*sin(theta);
B_x=posx+(auxdata.l_axes+auxdata.l_front).*cos(theta)+auxdata.b_width.*sin(theta);
C_x=posx-auxdata.l_rear.*cos(theta)+auxdata.b_width.*sin(theta);
D_x=posx-auxdata.l_rear.*cos(theta)-auxdata.b_width.*sin(theta);

f_p_A_x=(-tanh(50*(A_x-0.1))+tanh(50*(A_x-auxdata.SL+0.1)))/2.*auxdata.SW;
f_p_B_x=(-tanh(50*(B_x-0.1))+tanh(50*(B_x-auxdata.SL+0.1)))/2.*auxdata.SW;
f_p_C_x=(-tanh(50*(C_x-0.1))+tanh(50*(C_x-auxdata.SL+0.1)))/2.*auxdata.SW;
f_p_D_x=(-tanh(50*(D_x-0.1))+tanh(50*(D_x-auxdata.SL+0.1)))/2.*auxdata.SW;

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

c=[curvature_dot, auxdata.CL-A_y, auxdata.CL-B_y, auxdata.CL-C_y, auxdata.CL-D_y, A_y-f_p_A_x, B_y-f_p_B_x, C_y-f_p_C_x, D_y-f_p_D_x, AreaO-AreaRef, AreaE-AreaRef ];
%------------- END OF CODE --------------


function cr=avrc_unscaled(x,u,p,t,data)
% avrc_unscaled - Returns the rate constraint algebraic function where [xrl url] =<
% avrc(x,u,p,t) =< [xru uru]
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  cr=avrc_unscaled(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%   data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    cr - constraint function
%
%
%------------- BEGIN CODE --------------


[ cr ] = addRateConstraint( x,u,p,t,data );

%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
%------------- BEGIN CODE --------------
varargin=varargin{1};

auxdata = vdat.auxdata;

posyf = xf(2);
thetaf = xf(4);

A_y=posyf+(auxdata.l_axes+auxdata.l_front).*sin(thetaf)+auxdata.b_width.*cos(thetaf);
B_y=posyf+(auxdata.l_axes+auxdata.l_front).*sin(thetaf)-auxdata.b_width.*cos(thetaf);
C_y=posyf-auxdata.l_rear.*sin(thetaf)-auxdata.b_width.*cos(thetaf);
D_y=posyf-auxdata.l_rear.*sin(thetaf)+auxdata.b_width.*cos(thetaf);

bc=[A_y; B_y; C_y; D_y];

%------------- END OF CODE --------------
% When adpative time interval add constraint on time
%------------- BEGIN CODE --------------
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if ((strcmp(options.transcription,'hpLGR')) || (strcmp(options.transcription,'globalLGR')))  && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

function stageCost=L(x,xr,u,ur,p,t,vdat)

% L - Returns the stage cost.
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    if ~isempty(xr)
        xr=scale_variables_back( xr, vdat.Xscale, vdat.Xshift );
    end
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if ~isempty(ur)
        ur=scale_variables_back( ur, vdat.Uscale, vdat.Ushift );
    end
    stageCost=L_unscaled(x,u,p,t,vdat);
else
    stageCost=L_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,t0,tf,vdat) 

% E - Returns the boundary value cost
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift )';
    xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift )';
    u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift )';
    uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift )';
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
else
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
end


%------------- END OF CODE --------------


function dx = f(x,u,p,t,vdat)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    dx = f_unscaled(x,u,p,t,vdat);
    dx= scale_variables( dx, vdat.Xscale, 0 );
else
    dx = f_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
%     t=scale_variables_back( t, vdat.Tscale, vdat.Tshift );
    c = g_unscaled(x,u,p,t,vdat);
else
    c = g_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function cr=avrc(x,u,p,t,vdat)
% avrc - Returns the rate constraint algebraic function where [xrl url] =< avrc(x,u,p,t) =< [xru uru]
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    cr = avrc_unscaled(x,u,p,t,vdat);
else
    cr = avrc_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
% Warp function
%------------- BEGIN CODE --------------
bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
if isfield(vdat,'Xscale')
    if ~isempty(bc)
        x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift );
        xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift );
        u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift );
        uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift );
        if isfield(vdat,'Pscale')
            p=scale_variables_back( p', vdat.Pscale, vdat.Pshift );
        end
        bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
    end
end
%------------- END OF CODE --------------