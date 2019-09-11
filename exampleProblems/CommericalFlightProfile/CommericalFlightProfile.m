function [problem,guess] = CommericalFlightProfile
% CommericalFlightProfile - Minimum-Fuel Flight Profile for Commercial Aircraft
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
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
% Plant model name, used for Adigator
InternalDynamics=@CommericalFlightProfile_Dynamics_Internal;
SimDynamics=@CommericalFlightProfile_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_CommericalFlightProfile;

% Thrust model
x=[0 0.5 1];
y=[0 0.48 1];
FFModel=pchip(x,y);

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
auxdata.FFModel=FFModel;

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
auxdata.obs_npos_1=775000;
auxdata.obs_epos_1=842000;
auxdata.obs_r_1=40000;

auxdata.obs_npos_2=545000;
auxdata.obs_epos_2=442000;
auxdata.obs_r_2=30000;

auxdata.obs_npos_3=505000;
auxdata.obs_epos_3=842000;
auxdata.obs_r_3=20000;


auxdata.obs_npos_4=105000;
auxdata.obs_epos_4=72000;
auxdata.obs_r_4=60000;

auxdata.obs_npos_5=445000;
auxdata.obs_epos_5=0;
auxdata.obs_r_5=70000;

auxdata.obs_npos_6=600000;
auxdata.obs_epos_6=200000;
auxdata.obs_r_6=50000;

auxdata.obs_npos_7=100000;
auxdata.obs_epos_7=500000;
auxdata.obs_r_7=100000;

auxdata.obs_npos_8=300000;
auxdata.obs_epos_8=300000;
auxdata.obs_r_8=60000;

auxdata.obs_npos_8=800000;
auxdata.obs_epos_8=400000;
auxdata.obs_r_8=90000;

auxdata.obs_npos_9=300000;
auxdata.obs_epos_9=900000;
auxdata.obs_r_9=30000;


auxdata.obs_npos_10=600000;
auxdata.obs_epos_10=800000;
auxdata.obs_r_10=80000;
%%


%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=7000;     
problem.time.tf_max=10000; 
guess.tf=8500;

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
problem.states.xErrorTol_local=[1 1 1 0.5 deg2rad(5) deg2rad(5) 0.1*auxdata.g];
problem.states.xErrorTol_integral=[1 1 1 0.5 deg2rad(5) deg2rad(5) 0.1*auxdata.g];

% State constraint error bounds
problem.states.xConstraintTol=[1 1 1 0.5 deg2rad(5) deg2rad(5) 0.1*auxdata.g];
problem.states.xrConstraintTol=[1 1 1 0.5 deg2rad(5) deg2rad(5) 0.1*auxdata.g];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[Hf xf yf Vtasf gammaf chif W_min]; 
problem.states.xfu=[Hf xf yf Vtasf gammaf chif W0];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[H_max H_max];
guess.states(:,2)=[x0 xf];
guess.states(:,3)=[y0 yf];
guess.states(:,4)=[Vtas_max Vtas_max];
guess.states(:,5)=[gamma0 gammaf];
guess.states(:,6)=[chi0 chif];
guess.states(:,7)=[W0 (18000-2000)*auxdata.g];

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
problem.inputs.url=[-deg2rad(2) -deg2rad(10) -inf];
problem.inputs.uru=[deg2rad(2) deg2rad(10) inf];

% Input constraint error bounds
problem.inputs.uConstraintTol=[deg2rad(0.5) deg2rad(0.5) 0.1];
problem.inputs.urConstraintTol=[deg2rad(0.5) deg2rad(0.5) 0.1];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[alpha0  alpha0];
guess.inputs(:,2)=[phi_0  phi_0];
guess.inputs(:,3)=[1  1];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu

problem.constraints.ng_eq=0;
problem.constraints.gTol_eq=[];

problem.constraints.gl=[0 0 0 0 0 0 0 0 0 0];
problem.constraints.gu=[inf inf inf inf inf inf inf inf inf inf];
problem.constraints.gTol_neq=[(auxdata.obs_r_1+5)^2-(auxdata.obs_r_1)^2 (auxdata.obs_r_2+5)^2-(auxdata.obs_r_2)^2 (auxdata.obs_r_3+5)^2-(auxdata.obs_r_3)^2 (auxdata.obs_r_4+5)^2-(auxdata.obs_r_4)^2 (auxdata.obs_r_5+5)^2-(auxdata.obs_r_5)^2 (auxdata.obs_r_6+5)^2-(auxdata.obs_r_6)^2 (auxdata.obs_r_7+5)^2-(auxdata.obs_r_7)^2 (auxdata.obs_r_8+5)^2-(auxdata.obs_r_8)^2 (auxdata.obs_r_9+5)^2-(auxdata.obs_r_9)^2 (auxdata.obs_r_10+5)^2-(auxdata.obs_r_10)^2]/1e04;

% problem.constraints.gl=[0 0];
% problem.constraints.gu=[inf inf];
% problem.constraints.gTol_neq=[(auxdata.obs_r_1+5)^2-(auxdata.obs_r_1)^2 (auxdata.obs_r_4+5)^2-(auxdata.obs_r_4)^2]/1e04;


% % problem.constraints.g_neq_ActiveTime{1}=[guess.tf/2 guess.tf];
% problem.constraints.g_neq_ActiveTime{1}=[0 guess.tf];
% problem.constraints.g_neq_ActiveTime{2}=[];
% problem.constraints.g_neq_ActiveTime{3}=[];
% problem.constraints.g_neq_ActiveTime{4}=[0 guess.tf];
% % problem.constraints.g_neq_ActiveTime{4}=[0 guess.tf/2];
% problem.constraints.g_neq_ActiveTime{5}=[];



% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[0];
problem.constraints.bu=[0];
problem.constraints.bTol=[0.1];

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;

% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,vdat)

% L_unscaled - Returns the stage cost.
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

%------------- BEGIN CODE --------------


stageCost = 0*t;

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 

% E_unscaled - Returns the boundary value cost
%
% Syntax:  boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 
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

boundaryCost= tf;

%------------- END OF CODE --------------




function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
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
    if ((strcmp(options.discretization,'hpLGR')) || (strcmp(options.discretization,'globalLGR')))  && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leave the following unchanged! %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    if ~isempty(ur)
        ur=scale_variables_back( ur, vdat.Uscale, vdat.Ushift );
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost=0*t;
    else
        stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
    end
else
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost=0*t;
    else
        stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
    end
end

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,t0,tf,vdat) 

% E - Returns the boundary value cost
% Warp function
%------------- BEGIN CODE --------------
if isfield(vdat,'Xscale')
    x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift );
    xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift );
    u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift );
    uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p', vdat.Pscale, vdat.Pshift );
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        boundaryCost=sum(sum(p(:,end-vdat.mode.np*2+1:end)));
    else
        boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
    end
else
    if strcmp(vdat.mode.currentMode,'Feasibility')
        boundaryCost=sum(sum(p(:,end-vdat.mode.np*2+1:end)));
    else
        boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
    end
end


%------------- END OF CODE --------------


function dx = f(x,u,p,t,vdat)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% Warp function
%------------- BEGIN CODE --------------
f_unscaled=vdat.InternalDynamics;
if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    dx = f_unscaled(x,u,p,t,vdat);
    dx = scale_variables( dx, vdat.Xscale, 0 );
else
    dx = f_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------

function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
g_unscaled=vdat.InternalDynamics;
ng_group=nargout(g_unscaled);
if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
end

if ng_group==1
    c=[];
elseif ng_group==2
    [~,c] = g_unscaled(x,u,p,t,vdat);
else
    [~,ceq,cneq] = g_unscaled(x,u,p,t,vdat);
    c=[ceq cneq];
end

if isfield(vdat,'gFilter')
    c(:,vdat.gFilter)=[];
end

if strcmp(vdat.mode.currentMode,'Feasibility')
    c=[c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np) c+p(:,end-vdat.mode.np+1:end)];
end


function [dx,c] = fg(x,u,p,t,vdat)
% fg - Returns the ODE right hand side where x'= f(x,u,p,t) and the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
fg_unscaled=vdat.InternalDynamics;
ng_group=nargout(fg_unscaled);

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
end

if ng_group==1
    c=[];
elseif ng_group==2
    [dx,c]=fg_unscaled(x,u,p,t,vdat);
else
    [dx,ceq,cneq]=fg_unscaled(x,u,p,t,vdat);
    c=[ceq cneq];
end

if isfield(vdat,'Xscale')
    dx = scale_variables( dx, vdat.Xscale, 0 );
end

if isfield(vdat,'gFilter')
    c(:,vdat.gFilter)=[];
end

if strcmp(vdat.mode.currentMode,'Feasibility')
    c=[c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np) c+p(:,end-vdat.mode.np+1:end)];
end

%------------- END OF CODE --------------


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


%------------- END OF CODE ---------------------

function dx = f_unscaled(x,u,p,t,vdat)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% Warp function
%------------- BEGIN CODE --------------
Dynamics=vdat.InternalDynamics;
dx = Dynamics(x,u,p,t,vdat);

%------------- END OF CODE --------------

function c=g_unscaled(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
Dynamics=vdat.InternalDynamics;
ng_group=nargout(Dynamics);

if ng_group==1
    c=[];
elseif ng_group==2
    [~,c] = Dynamics(x,u,p,t,vdat);
else
    [~,ceq,cneq] = Dynamics(x,u,p,t,vdat);
    c=[ceq cneq];
end


%------------- END OF CODE ---------------------