function [problem,guess] = WindshearGoAround
%WindshearGoAround - Aircraft go around in the present of wind-shear problem
%
% The problem was originally presented by: 
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 1: Necessary conditions", Journal of Optimization Theory and Applications, 70(1), pp 1-23, 1991.
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 2: Multiple shooting and homotopy", Journal of Optimization Theory and Applications, 70(2), pp 223-254, 1991.
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% This implementation contains modifications
%
% Syntax:  [problem,guess] = WindshearGoAround
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
% 
%------------- BEGIN CODE --------------
% Aerodynamic coefficients and the fitting of look-up tables
a=6e-8;
b_var=-4e-11;
c=-log(25/30.6)*10^(-12);
e=6.280834899e-11;
d=-8.028808625e-8;
beta_0=0.3825;
beta_0_dot=0.2;
A_0=0.4456e05;
A_1=-0.2398e02;
A_2=0.1442e-01;
B_0=0.15523333333;
B_1=0.1236914764;
B_2=2.420265075;
C_0=0.7125;
C_1=6.087676573;
C_2=-9.027717451;
mg=150000;
g_var=3.2172e01;
delta=deg2rad(2);
rho=0.2203e-2;
S=0.1560e04;
alpha_star=deg2rad(12);
alpha_max=0.3002;

syms x
W1_1= double(sym2poly(-50+a*x^3+b_var*x^4));
W1_2= [zeros(1,3) double(sym2poly((x+500-2300)/40))];
W1_3= double(sym2poly(50-a*(4600-(x+4100))^3-b_var*(4600-(x+4100))^4));
W1_4= [0 0 0 0 50];
W1_poly = mkpp([0,500,4100,4600,20000],[W1_1;W1_2;W1_3;W1_4]);

W1_dot_1= double(sym2poly(diff(-50+a*x^3+b_var*x^4)));
W1_dot_2= [zeros(1,3) double(sym2poly(diff((x+500-2300)/40)))];
W1_dot_3= double(sym2poly(diff(50-a*(4600-(x+4100))^3-b_var*(4600-(x+4100))^4)));
W1_dot_4= [0 0 0 0];
W1_dot_poly = mkpp([0,500,4100,4600,20000],[W1_dot_1;W1_dot_2;W1_dot_3;W1_dot_4]);

xp=500:200:4100;
y=-51.*exp(-c.*(xp-2300).^4);
p = pchip(xp,y);
W2_1=double(sym2poly(d*x^3+e*x^4));
W2_2=[zeros(length(p.coefs),1) p.coefs];
W2_3=double(sym2poly(d*(4600-(x+4100))^3+e*(4600-(x+4100))^4));
W2_4=zeros(1,5);
W2_poly = mkpp([0,p.breaks,4600,20000],[W2_1;W2_2;W2_3;W2_4]);

dy(x)=diff(-51.*exp(-c.*(x-2300).^4));
dy=double(dy(xp));
p_dot = pchip(xp,dy);
W2_dot_1=double(sym2poly(diff(d*x^3+e*x^4)));
W2_dot_2=[p_dot.coefs];
W2_dot_3=double(sym2poly(diff(d*(4600-(x+4100))^3+e*(4600-(x+4100))^4)));
W2_dot_4=zeros(1,4);
W2_dot_poly = mkpp([0,p.breaks,4600,20000],[W2_dot_1;W2_dot_2;W2_dot_3;W2_dot_4]);

beta_poly=mkpp([0, (1-beta_0)/beta_0_dot, 1000],[beta_0_dot beta_0; 0 1]);

syms x
CL_1= double(sym2poly(C_0+C_1*(x-alpha_star)));
CL_2= double(sym2poly(C_0+C_1*(x+alpha_star)+C_2*((x+alpha_star)-alpha_star)^2));
C_L_poly= mkpp([-alpha_star, alpha_star, alpha_max ],[[0 CL_1];CL_2]);


% Store relevant data
auxdata.a=a;
auxdata.b=b_var;
auxdata.c=c;
auxdata.e=e;
auxdata.d=d;
auxdata.beta_0=beta_0;
auxdata.beta_0_dot=beta_0_dot;
auxdata.A_0=A_0;
auxdata.A_1=A_1;
auxdata.A_2=A_2;
auxdata.B_0=B_0;
auxdata.B_1=B_1;
auxdata.B_2=B_2;
auxdata.C_0=C_0;
auxdata.C_1=C_1;
auxdata.C_2=C_2;
auxdata.g=g_var;
auxdata.m=mg/g_var;
auxdata.delta=delta;
auxdata.rho=rho;
auxdata.S=S;
auxdata.alpha_star=alpha_star;
auxdata.alpha_max=alpha_max;

auxdata.W1_poly=W1_poly;
auxdata.W2_poly=W2_poly;
auxdata.W1_dot_poly=W1_dot_poly;
auxdata.W2_dot_poly=W2_dot_poly;
auxdata.beta_poly=beta_poly;
auxdata.C_L_poly=C_L_poly;

% initial conditions
pos0 = 0;
alt0 = 600; % Initial altitude (ft)
speed0 = 239.7; % Initial speed (ft/s)
fpa0 = -0.03925; % Initial flight path angle (rad)
fpaf = 0.12969; % Final flight path angle (rad)
alpha0 = 0.1283; % Initial angle of attack (rad)

% variable bounds
posmin = 0; posmax = 1e04;
altmin = 1e02; altmax = 1000;
speedmin = 0; speedmax = inf;
fpamin = -inf; fpamax = inf;
alphamin = -alpha_max; alphamax = alpha_max;
umin = -0.05236; umax = 0.05236;

%%

%------------- BEGIN CODE --------------
% Plant model name
InternalDynamics=@WindshearGoAround_Dynamics_Internal;
SimDynamics=@WindshearGoAround_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_WindshearGoAround;

%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=40;     
problem.time.tf_max=40; 
guess.tf=40;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=300;
problem.parameters.pu=inf;
guess.parameters=502;

% Initial conditions for system.
problem.states.x0=[pos0 alt0 speed0 fpa0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[pos0 alt0 speed0 fpa0]; 
problem.states.x0u=[pos0 alt0 speed0 fpa0]; 

% State bounds. xl=< x <=xu
problem.states.xl=[posmin altmin speedmin fpamin]; 
problem.states.xu=[posmax altmax speedmax fpamax]; 

% State rate bounds. xrl=< x <=xru
problem.states.xrl=[-inf -inf -inf -inf]; 
problem.states.xru=[inf inf inf inf]; 

% State error bounds
problem.states.xErrorTol_local=[1 0.5 0.1 deg2rad(0.5)];
problem.states.xErrorTol_integral=[1 0.5 0.1 deg2rad(0.5)];

% State constraint error bounds
problem.states.xConstraintTol=[1 0.5 0.1 deg2rad(0.5)];
problem.states.xrConstraintTol=[1 0.5 0.1 deg2rad(0.5)];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[posmin altmin speedmin fpaf]; 
problem.states.xfu=[posmax altmax speedmax fpaf];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[pos0 9000];
guess.states(:,2)=[alt0 850];
guess.states(:,3)=[speed0 speed0];
guess.states(:,4)=[fpa0 fpaf];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[alphamin];
problem.inputs.uu=[alphamax];

problem.inputs.u0l=[alphamin];
problem.inputs.u0u=[alphamax];

% Input rate bounds
problem.inputs.url=[umin];
problem.inputs.uru=[umax];

% Input constraint error bounds
problem.inputs.uConstraintTol=[deg2rad(0.5)];
problem.inputs.urConstraintTol=[deg2rad(0.5)];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[alpha0 alpha0];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.ng_eq=0;
problem.constraints.gTol_eq=[];

problem.constraints.gl=[0];
problem.constraints.gu=[inf];
problem.constraints.gTol_neq=[1e-2];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];
problem.constraints.bTol=[];

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

boundaryCost=-p(1);

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
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

