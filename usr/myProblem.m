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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
%
%Initial Time. t0<tf
problem.time.t0=0;

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

% State rate bounds. xrl=< x <=xru
problem.states.xrl=[x1dot_lowerbound ... xndot_lowerbound]; 
problem.states.xru=[x1dot_upperbound ... xndot_upperbound]; 

% State error bounds
problem.states.xErrorTol=[eps_x1 ... eps_xn]; 

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

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[g1_lowerbound g2_lowerbound ...];
problem.constraints.gu=[g1_upperbound g2_upperbound ...];

% Path constraint error bounds
problem.constraints.gTol=[eps_g1_bounds eps_g2_bounds ...];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[b1_lowerbound b2_lowerbound ...];
problem.constraints.bu=[b1_upperbound b2_upperbound ...];

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;

% Plant model name, for use with Adigator
problem.data.plantmodel = '...';

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

%Define states and setpoints
x1 = x(:,1); xr1=xr(:,1);
%...
xn=x(:,n); xrn=xr(:,n);

%Define inputs
u1 = u(:,1);ur1=ur(:,1);
% ...
um = u(:,m);urm=ur(:,m);

stageCost = ...;

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

boundaryCost=...;

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

%Stored data
auxdata = data.auxdata;

%Define states
x1 = x(:,1);
%...
xn=x(:,n),

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);


%Define ODE right-hand side

dx(:,1) = f1(x1,..xn,u1,..um,p,t);
%...
dx(:,n) = fn(x1,..xn,u1,..um,p,t);


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

%Define states
x1 = x(:,1);
%...
xn=x(:,n),

%Define inputs
u1 = u(:,1);
% ...
um = u(:,m);

c(:,1)=g1(x1,...,u1,...p,t);
c(:,2)=g2(x1,...,u1,...p,t);

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

% when having rate constraints (i.e. not all equal to +/-inf)
[ cr ] = addRateConstraint( x,u,p,t,data );

% when not having rate constraints
% [ cr ] = [];

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