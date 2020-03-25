function [problem,guess] = ShuttleReentryTrajectory
%ShuttleReentryTrajectory - Space Shuttle Re-entry Trajectory Problem
%
% The problem was adapted from Example 6.1 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
%
% Syntax:  [problem,guess] = MinTimeClimbBryson
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


% Boundary Conditions
h_0 = 260000;
phi_0 = 0;
theta_0 = 0;
v_0 = 25600;
gamma_0 = -1*pi/180;
psi_0 = 90*pi/180;

h_f = 80000;
v_f = 2500;
gamma_f = -5*pi/180;

% variable simple bounds
h_min = 0; h_max = h_0;
phi_min = -pi; phi_max = pi;
theta_min = -89*pi/180; theta_max = 89*pi/180;
v_min = 1; v_max = v_0*1.1;
gamma_min = -89*pi/180; gamma_max = 89*pi/180;
psi_min = -90*pi/180; psi_max = 90*pi/180;
alpha_min = -90*pi/180; alpha_max = 90*pi/180;
beta_min = -89*pi/180; beta_max = 1*pi/180;

%%

% Plant model name, used for Adigator
InternalDynamics=@ShuttleReentryTrajectory_Dynamics_Internal;
SimDynamics=@ShuttleReentryTrajectory_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_ShuttleReentryTrajectory;

%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=0;     
problem.time.tf_max=3000; 
guess.tf=2000;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[h_0 phi_0 theta_0 v_0 gamma_0 psi_0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[h_0 phi_0 theta_0 v_0 gamma_0 psi_0]; 
problem.states.x0u=[h_0 phi_0 theta_0 v_0 gamma_0 psi_0]; 

% State bounds. xl=< x <=xu
problem.states.xl=[h_min phi_min theta_min v_min gamma_min psi_min];
problem.states.xu=[h_max phi_max theta_max v_max gamma_max psi_max];


% State error bounds
problem.states.xErrorTol_local=[1 deg2rad(0.5) deg2rad(0.5) 0.1 deg2rad(0.5) deg2rad(0.5)];
problem.states.xErrorTol_integral=[10 deg2rad(0.5) deg2rad(0.5) 0.1 deg2rad(0.5) deg2rad(0.5)];

% State constraint error bounds
problem.states.xConstraintTol=[1 deg2rad(0.5) deg2rad(0.5) 0.1 deg2rad(0.5) deg2rad(0.5)];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[h_f phi_min theta_min v_f gamma_f psi_min]; 
problem.states.xfu=[h_f phi_max theta_max v_f gamma_f psi_max];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[h_0 h_f];
guess.states(:,2)=[phi_0 phi_0+10*pi/180];
guess.states(:,3)=[theta_0 theta_0+10*pi/180];
guess.states(:,4)=[v_0 v_f];
guess.states(:,5)=[gamma_0 gamma_f];
guess.states(:,6)=[psi_0 (psi_max+psi_min)/2];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[alpha_min beta_min];
problem.inputs.uu=[alpha_max beta_max];

problem.inputs.u0l=[alpha_min beta_min];
problem.inputs.u0u=[alpha_max beta_max];


% Input constraint error bounds
problem.inputs.uConstraintTol=[deg2rad(0.5) deg2rad(0.5)];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[alpha_min alpha_max];
guess.inputs(:,2)=[beta_min beta_max];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.ng_eq=0;
problem.constraints.gTol_eq=[];

problem.constraints.gl=[];
problem.constraints.gu=[];
problem.constraints.gTol_neq=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];
problem.constraints.bTol=[];


% store the necessary problem parameters used in the functions
problem.data.a0 = -0.20704;
problem.data.a1 = 0.029244;
problem.data.mu = 0.14076539e17;
problem.data.Re = 20902900;
problem.data.S  = 2690;
problem.data.b0 = 0.07854;
problem.data.b1 = -0.0061592;
problem.data.b2 = 0.000621408;
problem.data.hr = 23800;
problem.data.rho0 = 0.002378;
problem.data.mass = 203000/32.174;


% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.data.functions_unscaled=problem.functions_unscaled;
problem.data.ng_eq=problem.constraints.ng_eq;
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

boundaryCost=-xf(3);

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
bc=[];
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

