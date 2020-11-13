function [problem,guess] = myProblem
%myProblem - Template file for optimal control problem definition
%
%Syntax:  [problem,guess] = myProblem
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
% Plant model name, provide in the format of function handle
InternalDynamics=@myProblem_Dynamics_Internal; 
SimDynamics=@myProblem_Dynamics_Sim;

% Analytic derivative files (optional), provide in the format of function handle
problem.analyticDeriv.gradCost=@gradCost_myProblem;
problem.analyticDeriv.hessianLagrangian=@hessianLagrangian_myProblem;
problem.analyticDeriv.jacConst=@jacConst_myProblem;

% (optional) customized call back function
%problem.callback=@callback_myProblem;

% Settings file
problem.settings=@settings_myProblem;

%Initial Time. t0<tf
problem.time.t0_min=initial_time_min;
problem.time.t0_max=final_time_max;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=final_time_min;     
problem.time.tf_max=final_time_max; 
guess.tf=final_time_guess;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[p1_lowerbound ...];
problem.parameters.pu=[p1_upperbound ...];
guess.parameters=[p1_guess p2_guess ...];

% Initial conditions for system.
problem.states.x0=[x1(t0) ... xn(t0)];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[x1(t0)_lowerbound ... xn(t0)_lowerbound]; 
problem.states.x0u=[x1(t0)_upperbound ... xn(t0)_upperbound]; 

% State bounds. xl=< x <=xu
problem.states.xl=[x1_lowerbound ... xn_lowerbound];
problem.states.xu=[x1_upperbound ... xn_upperbound];

% State rate bounds. xrl=< x_dot <=xru
problem.states.xrl=[x1dot_lowerbound ... xndot_lowerbound]; 
problem.states.xru=[x1dot_upperbound ... xndot_upperbound]; 

% State error bounds
problem.states.xErrorTol_local=[eps_loc_x1 ... eps_loc_xn]; 
problem.states.xErrorTol_integral=[eps_int_x1 ... eps_int_xn]; 

% State constraint error bounds
problem.states.xConstraintTol=[eps_x1_bounds ... eps_xn_bounds];
problem.states.xrConstraintTol=[eps_x1dot_bounds ... eps_xndot_bounds];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[x1(tf)_lowerbound ... xn(tf)_lowerbound]; 
problem.states.xfu=[x1(tf)_upperbound ... xn(tf)_upperbound];

% Guess the state trajectories with [x0 ... xf]
guess.time=[t0 ... tf];
guess.states(:,1)=[x1(t0) ... x1(tf)];
% ...
guess.states(:,n)=[xn(t0) ... xn(tf)];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[u1_lowerbound ... um_lowerbound];
problem.inputs.uu=[u1_upperbound ... um_upperbound];

% Bounds on the first control action
problem.inputs.u0l=[u1_lowerbound(t0) ... um_lowerbound(t0)];
problem.inputs.u0u=[u1_upperbound(t0) ... um_upperbound(t0)];

% Input rate bounds
problem.inputs.url=[u1dot_lowerbound ... umdot_lowerbound]; 
problem.inputs.uru=[u1dot_upperbound ... umdot_upperbound]; 

% Input constraint error bounds
problem.inputs.uConstraintTol=[eps_u1_bounds ... eps_um_bounds];
problem.inputs.urConstraintTol=[eps_u1dot_bounds ... eps_umdot_bounds];

% Guess the input sequences with [u0 ... uf]
guess.inputs(:,1)=[u1(t0) ... u1(tf)]];
%...
guess.inputs(:,m)=[um(t0) ... um(tf)]];

% Path constraint function 
problem.constraints.ng_eq=0; % number of quality constraints in format of g(x,u,p,t) == 0
problem.constraints.gTol_eq=[eps_g1_bounds eps_g2_bounds ...]; % equality cosntraint error bounds

problem.constraints.gl=[g1_lowerbound g2_lowerbound ...]; % Lower ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gu=[g1_upperbound g2_upperbound ...]; % Upper ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gTol_neq=[eps_g1_bounds eps_g2_bounds ...]; % inequality constraint error bounds

% OPTIONAL: define the time duration each constraint will be active, for
% example (for ECH enabled in setings)
problem.constraints.gActiveTime{1}=[guess.tf/2 guess.tf];
problem.constraints.gActiveTime{2}=[];
...
problem.constraints.gActiveTime{5}=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[b1_lowerbound b2_lowerbound ...];
problem.constraints.bu=[b1_upperbound b2_upperbound ...];
problem.constraints.bTol=[eps_b1_bounds eps_b2_bounds ...]; 

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;
% optional setting for automatic regularization
problem.data.penalty.values=[weight_1, weight_2, ... weight_n];
problem.data.penalty.i=1; %starting weight



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

function stageCost=L_unscaled(x,xr,u,ur,p,t,data)

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

%Define states and setpoints
x1 = x(:,1); 
%...
xn=x(:,n); 

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);

stageCost = ...;

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

boundaryCost=...;

%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,data,varargin)

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
% Leave it here
varargin=varargin{1};
%------------- BEGIN CODE --------------
bc(1,:)=b1(x0,xf,u0,uf,p,tf);
bc(2,:)=b2(x0,xf,u0,uf,p,tf);
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
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

