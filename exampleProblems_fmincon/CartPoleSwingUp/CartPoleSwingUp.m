function [problem,guess] = CartPoleSwingUp
%CartPoleSwingUp - Cart Pole Swing-up Problem
%
% The problem was adapted from Section 6 of
% Kelly M. An introduction to trajectory optimization: How to do your own direct collocation. SIAM Review. 2017 Nov 6;59(4):849-904.
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
InternalDynamics=@CartPoleSwingUp_Dynamics_Internal;
SimDynamics=@CartPoleSwingUp_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=@gradCost_CartPoleSwingUp;
problem.analyticDeriv.hessianLagrangian=@hessianLagrangian_CartPoleSwingUp;
problem.analyticDeriv.jacConst=@jacConst_CartPoleSwingUp;

% Settings file
problem.settings=@settings_CartPoleSwingUp;

% Initial time. t0<tf
problem.time.t0_min=0;     
problem.time.t0_max=0; 
guess.t0=0;


% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=2;     
problem.time.tf_max=2; 
guess.tf=2;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system
problem.states.x0=[0 0 0 0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[0 0 0 0]; 
problem.states.x0u=[0 0 0 0]; 

% State bounds. xl=< x <=xu
problem.states.xl=[-inf -inf -2 -2*pi];
problem.states.xu=[inf inf 2 2*pi];

% State error bounds
problem.states.xErrorTol_local=[1e-5 1e-5 1e-5 1e-5];
problem.states.xErrorTol_integral=[1e-5 1e-5 1e-5 1e-5];

% State constraint error bounds
problem.states.xConstraintTol=[1e-5 1e-5 1e-5 1e-5];


% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[0 0 1 pi]; 
problem.states.xfu=[0 0 1 pi];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[0 0];
guess.states(:,2)=[0 0];
guess.states(:,3)=[0 1];
guess.states(:,4)=[0 pi];



% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;

% Input bounds
problem.inputs.ul=[-20];
problem.inputs.uu=[20];

% Bounds on the first control action
problem.inputs.u0l=[-20];
problem.inputs.u0u=[20];


% Input constraint error bounds
problem.inputs.uConstraintTol=[0.1];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[0, 0];   


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
problem.data.g = 9.81;
problem.data.m1 = 1;
problem.data.m2 = 0.3;
problem.data.L  = 0.5;

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


stageCost = u(:,1).*u(:,1);

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat) 

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

boundaryCost=0;

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
            bc=[bc;diff(t_segment)'];
        end
    end
end

%------------- END OF CODE --------------

