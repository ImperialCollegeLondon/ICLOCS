function [problem,guess] = LowThrustOrbitTransfer
%LowThrustOrbitTransfer - Low Thrust Orbit Transfer Problem
%
% The problem was adapted from Example 6.3 from
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
%
% Syntax:  [problem,guess] = LowThrustOrbitTransfer
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

% initial conditions
p_0 = 21837080.052835; 
f_0 = 0;
g_0 = 0;
h_0 = -0.25396764647494;
k_0 = 0;
L_0 = pi;
w_0 = 1;

% terminal conditions
p_f = 40007346.015232;

% variable lower bounds
p_min = 20000000; 
f_min = -1; 
g_min = -1; 
h_min = -1; 
k_min = -1; 
L_min = L_0; 
w_min = 0.1; 
u_r_min = -1; 
u_t_min = -1; 
u_h_min = -1; 
tau_min = -50; 

% variable upper bounds
p_max = 60000000;
f_max = +1;
g_max = +1;
h_max = +1;
k_max = +1;
L_max = 10*2*pi;
w_max = w_0;
u_r_max = +1;
u_t_max = +1;
u_h_max = +1;
tau_max = 0;

% boundary conditions
bc1=0.73550320568829;
bc2=0.61761258786099;

%%

% Plant model name, used for Adigator
InternalDynamics=@LTOT_Dynamics_Internal;
SimDynamics=@LTOT_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_LowThrustOrbitTransfer;

%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=50000;     
problem.time.tf_max=100000; 
guess.tf=90000;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=tau_min;
problem.parameters.pu=tau_max;
guess.parameters=-8;

% Initial conditions for system.
problem.states.x0=[p_0 f_0 g_0 h_0 k_0 L_0 w_0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[p_0 f_0 g_0 h_0 k_0 L_0 w_0];
problem.states.x0u=[p_0 f_0 g_0 h_0 k_0 L_0 w_0];

% State bounds. xl=< x <=xu
problem.states.xl=[p_min f_min g_min h_min k_min L_min w_min];
problem.states.xu=[p_max f_max g_max h_max k_max L_max w_max];

% State error bounds
problem.states.xErrorTol_local=[100 1 1 1 1 1 1];
problem.states.xErrorTol_integral=[100 1 1 1 1 1 1];

% State constraint error bounds
problem.states.xConstraintTol=[1 1 1 1 1 1 1];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[p_f f_min g_min h_min k_min L_min w_min];
problem.states.xfu=[p_f f_max g_max h_max k_max L_max w_max];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[u_r_min u_t_min u_h_min];
problem.inputs.uu=[u_r_max u_t_max u_h_max];

problem.inputs.u0l=[u_r_min u_t_min u_h_min];
problem.inputs.u0u=[u_r_max u_t_max u_h_max];

% Input constraint error bounds
problem.inputs.uConstraintTol=[1 1 1];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.ng_eq=1;
problem.constraints.gTol_eq=[1];

problem.constraints.gl=[];
problem.constraints.gu=[];
problem.constraints.gTol_neq=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[bc1.^2 bc2.^2 0 -3];
problem.constraints.bu=[bc1.^2 bc2.^2 0 0];
problem.constraints.bTol=[1 1 1 1];


% store the necessary problem parameters used in the functions
problem.data.g0 = 32.174;
problem.data.Isp = 450; 
problem.data.T = 4.446618e-03; 
problem.data.mu = 1.407645794e16; 
problem.data.Re = 20925662.73; 
problem.data.J2 = 1082.639e-06;
problem.data.J3 = -2.565e-06;
problem.data.J4 = -1.608e-06;

% Obtain guess of states and input sequences with ode solve
[guess.time,guess.states] = ode45(@(t,x) odeInitialGuess(t,x,SimDynamics,problem.data), linspace(guess.t0,guess.tf,1000), problem.states.x0);
guess.inputs(:,1)=-sind(0)*ones(size(guess.time));
guess.inputs(:,2)=cosd(0)*cosd(30)*ones(size(guess.time));
guess.inputs(:,3)=-cosd(0)*sind(30)*ones(size(guess.time));

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

boundaryCost=-xf(7);

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
bc=[xf(2)^2+xf(3)^2;xf(4)^2+xf(5)^2;xf(2)*xf(4)+xf(3)*xf(5);xf(3)*xf(4)-xf(5)*xf(2)];
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

