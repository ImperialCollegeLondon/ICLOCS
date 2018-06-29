function [problem,guess] = F50MinTimeFlight
% F50MinTimeFlight - Minimum Time Flight Profile for Commercial Aircraft
%
% The aerodynamic and propulsion data are obtained from
% "Performance model Fokker 50", Delft University of Technology, 2010
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
problem.data.plantmodel = 'F50MinTimeFlightPlant';

% Constants
auxdata.g=9.81;
auxdata.ktomps=0.514444;
auxdata.ftom=0.3048;
auxdata.R=287.15;
auxdata.ps=101325;
auxdata.rhos=1.225;
auxdata.Ts=288.15;
auxdata.dTdH=-0.0065;
auxdata.nprop=1020;
auxdata.Dprop=3.66;
auxdata.kappa=1.4;
auxdata.alpha0=-2.32*pi/180;
auxdata.clalpha=0.095;
auxdata.k_cd=0.033;
auxdata.cd0=0.0233;
auxdata.S=70;
auxdata.alphat=0*pi/180;

% Look up tables
auxdata.a1p=1.5456*10^-9;
auxdata.a2p=-3.1176*10^-7;
auxdata.a3p=1.9477*10^-5;
auxdata.b1p=-2.0930*10^-5;
auxdata.b2p=4.1670*10^-3;
auxdata.b3p=-0.40739;
auxdata.c1p=0.088835;
auxdata.c2p=-12.855;
auxdata.c3p=3077.7;
auxdata.Pidle=75;

% initial conditions
W0=18000*auxdata.g;
H0=1600*auxdata.ftom;
T0=auxdata.Ts+H0*auxdata.dTdH;
p0=auxdata.ps*(T0/auxdata.Ts)^(-auxdata.g/auxdata.dTdH/auxdata.R);
rho0=auxdata.rhos*(T0/auxdata.Ts)^(-(auxdata.g/auxdata.dTdH/auxdata.R+1));
Vtas0=CAS2TAS( auxdata.kappa, p0, rho0, auxdata.ps, auxdata.rhos, 190*auxdata.ktomps );
x0=0;y0=0;
alpha0=5*pi/180;
chi0=0;
gamma0=8*pi/180;

% variable bounds
W_min=(18000-4000)*auxdata.g;
H_min=H0;H_max=25000*auxdata.ftom;
T_Hmax=auxdata.Ts+H_max*auxdata.dTdH;
p_Hmax=auxdata.ps*(T_Hmax/auxdata.Ts)^(-auxdata.g/auxdata.dTdH/auxdata.R);
rho_Hmax=auxdata.rhos*(T_Hmax/auxdata.Ts)^(-(auxdata.g/auxdata.dTdH/auxdata.R+1));
Vtas_max=CAS2TAS( auxdata.kappa, p_Hmax, rho_Hmax, auxdata.ps, auxdata.rhos, 227*auxdata.ktomps );
Vtas_min=Vtas0;
alpha_max=10*pi/180;
gamma_max=60*pi/180;
throttle_max=1;throttle_min=0;
chi_max=pi;
chi_min=-pi;
phi_0=0;
phi_max=45*pi/180;
x_max=1200000;
y_max=1200000;

% terminal conditions
Hf=2000*auxdata.ftom;
T_Hf=auxdata.Ts+Hf*auxdata.dTdH;
p_Hf=auxdata.ps*(T_Hf/auxdata.Ts)^(-auxdata.g/auxdata.dTdH/auxdata.R);
rho_Hf=auxdata.rhos*(T_Hf/auxdata.Ts)^(-(auxdata.g/auxdata.dTdH/auxdata.R+1));
Vtasf=CAS2TAS( auxdata.kappa, p_Hf, rho_Hf, auxdata.ps, auxdata.rhos, 190*auxdata.ktomps );
xf=800000;
yf=900000;
chif=-3/4*pi;
gammaf=-3*pi/180;

% Parameter for the obstacle
auxdata.obs_npos=745000;
auxdata.obs_epos=842000;
auxdata.obs_r=30000;

%%
%Initial Time. t0<tf
problem.time.t0=0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=0;     
problem.time.tf_max=15000; 
guess.tf=10000;
% guess.tf=presol.T(end);

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[H0 x0 y0 Vtas0 gamma0 chi0 W0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[H0 x0 y0 Vtas0 gamma0 chi0 W0];
problem.states.x0u=[H0 x0 y0 Vtas0 gamma0 chi0 W0];

% State bounds. xl=< x <=xu
problem.states.xl=[H_min x0 y0 Vtas_min -gamma_max chi_min W_min]; 
problem.states.xu=[H_max x_max y_max Vtas_max gamma_max chi_max W0]; 

% State rate bounds. xrl=< x <=xru
problem.states.xrl=[-inf -inf -inf -inf -inf -inf -inf]; 
problem.states.xru=[inf inf inf inf inf inf inf]; 

% State error bounds
problem.states.xErrorTol=[5 1 1 0.5 deg2rad(0.5) deg2rad(0.5) 1];

% State constraint error bounds
problem.states.xConstraintTol=[5 1 1 0.5 deg2rad(0.5) deg2rad(0.5) 1];
problem.states.xrConstraintTol=[5 1 1 0.5 deg2rad(0.5) deg2rad(0.5) 1];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[Hf xf yf Vtasf gammaf chif W_min]; 
problem.states.xfu=[Hf xf yf Vtasf gammaf chif W0];

% Guess the state trajectories with [x0 xf]
guess.time=[0 guess.tf/2 guess.tf];
guess.states(:,1)=[H_max H_max H_max];
guess.states(:,2)=[x0 xf/2 xf];
guess.states(:,3)=[y0 yf/2 yf];
guess.states(:,4)=[Vtas_max Vtas_max Vtas_max];
guess.states(:,5)=[gamma0 gammaf gammaf];
guess.states(:,6)=[chi0 chi0 chif];
guess.states(:,7)=[W0 (18000-2000)*auxdata.g W_min];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.

problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[-alpha_max -phi_max throttle_min];
problem.inputs.uu=[alpha_max phi_max throttle_max];

problem.inputs.u0l=[-alpha_max -phi_max throttle_min];
problem.inputs.u0u=[alpha_max phi_max throttle_max];

% Input rate bounds
problem.inputs.url=[-deg2rad(2) -deg2rad(10) -0.2];
problem.inputs.uru=[deg2rad(2) deg2rad(10) 0.2];

% Input constraint error bounds
problem.inputs.uConstraintTol=[deg2rad(0.5) deg2rad(0.5) 0.1];
problem.inputs.urConstraintTol=[deg2rad(0.5) deg2rad(0.5) 0.1];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[alpha0 alpha0 alpha0];
guess.inputs(:,2)=[phi_0 phi_0 phi_0];
guess.inputs(:,3)=[1 1 1];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[0];
problem.constraints.gu=[inf];
problem.constraints.gTol=[0.1];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[0];
problem.constraints.bu=[0];

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc_unscaled,@b_unscaled};
problem.constraintErrorTol=[problem.constraints.gTol,problem.constraints.gTol,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,data)


% L - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
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

% E - Returns the boundary value cost
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
%
boundaryCost=tf;

%------------- END OF CODE --------------


function dx = f_unscaled(x,u,p,t,data)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
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
W_dot= calcWeight(V_cas,H,throttle);

dx = [H_dot, x_dot, y_dot, TAS_dot, gamma_dot, chi_dot, W_dot];


%------------- END OF CODE --------------


function c=g_unscaled(x,u,p,t,data)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
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

npos = x(:,2);
epos = x(:,3);

c=[(npos-data.auxdata.obs_npos).^2+(epos-data.auxdata.obs_epos).^2-data.auxdata.obs_r.^2];

%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
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
bc=[uf(2)];
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
%------------- BEGIN CODE --------------
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if strcmp(options.transcription,'hpLGR') && options.adaptseg==1 
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
%     xr=scale_variables( xr, 1./vdat.Xscale );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
%     ur=scale_variables( ur, 1./vdat.Xscale );
    stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
else
    stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
end

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,t0,tf,vdat) 

% E - Returns the boundary value cost
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
if isfield(vdat,'Xscale')
    x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift );
    xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift );
    u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift );
    uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p', vdat.Pscale, vdat.Pshift );
    end
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
else
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
end


%------------- END OF CODE --------------


function dx = f(x,u,p,t,vdat)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
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

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    dx = f_unscaled(x,u,p,t,vdat);
    dx= scale_variables( dx, vdat.Xscale, 0 );
else
    dx = f_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function cr=avrc(x,u,p,t,data)

% avrc - Returns the rate constraint algebraic function where [xrl url] =<
% avrc(x,u,p,t) =< [xru uru]
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  cr=avrc(x,u,p,t,data)
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

function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
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
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    c = g_unscaled(x,u,p,t,vdat);
else
    c = g_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
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