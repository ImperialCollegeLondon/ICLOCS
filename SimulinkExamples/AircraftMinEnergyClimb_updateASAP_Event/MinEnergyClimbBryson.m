function [problem,guess] = MinEnergyClimbBryson
%MinFuelClimbBryson - Supersonic Aircraft Minimum Fuel Climb Problem
%
% The problem was adapted from the supersonic aircraft minimum time-to-climb problem originally presented by
% A. E. Bryson, M. N. Desai, and W. C. Hoffman, "Energy-State Approximation in Performance Optimization of Supersonic Aircraft," Journal of Aircraft, Vol. 6, No. 6, November-December, 1969, pp. 481-488. 
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% With aerodynamic data and modifications to thrust data by:
% M.A. Patterson and A.V. Rao, "GPOPS-II: A General Purpose MATLAB Software for Solving Multiple-Phase Optimal Control Problems, User's Manual, 2.3 edition", 2016
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
% Plant model name, used for Adigator
problem.data.plantmodel = 'MinEnergyClimbBrysonPlant';

% Lookup Table 1: U.S. 1976 Standard Atmosphere (altitude, density and pressure)
AtomsData = [-2000 1.478e+00 3.479e+02
0 1.225e+00 3.403e+02
2000 1.007e+00 3.325e+02
4000 8.193e-01 3.246e+02
6000 6.601e-01 3.165e+02
8000 5.258e-01 3.081e+02
10000 4.135e-01 2.995e+02
12000 3.119e-01 2.951e+02
14000 2.279e-01 2.951e+02
16000 1.665e-01 2.951e+02
18000 1.216e-01 2.951e+02
20000 8.891e-02 2.951e+02
22000 6.451e-02 2.964e+02
24000 4.694e-02 2.977e+02
26000 3.426e-02 2.991e+02
28000 2.508e-02 3.004e+02
30000 1.841e-02 3.017e+02
32000 1.355e-02 3.030e+02
34000 9.887e-03 3.065e+02
36000 7.257e-03 3.101e+02
38000 5.366e-03 3.137e+02
40000 3.995e-03 3.172e+02
42000 2.995e-03 3.207e+02
44000 2.259e-03 3.241e+02
46000 1.714e-03 3.275e+02
48000 1.317e-03 3.298e+02
50000 1.027e-03 3.298e+02
52000 8.055e-04 3.288e+02
54000 6.389e-04 3.254e+02
56000 5.044e-04 3.220e+02
58000 3.962e-04 3.186e+02
60000 3.096e-04 3.151e+02
62000 2.407e-04 3.115e+02
64000 1.860e-04 3.080e+02
66000 1.429e-04 3.044e+02
68000 1.091e-04 3.007e+02
70000 8.281e-05 2.971e+02
72000 6.236e-05 2.934e+02
74000 4.637e-05 2.907e+02
76000 3.430e-05 2.880e+02
78000 2.523e-05 2.853e+02
80000 1.845e-05 2.825e+02
82000 1.341e-05 2.797e+02
84000 9.690e-06 2.769e+02
86000 6.955e-06 2.741e+02];

% Lookup Table 2: Aerodynamic Data 
MachLTAero = [0 0.4 0.8 0.9 1.0 1.2 1.4 1.6 1.8];
ClalphaLT = [3.44 3.44 3.44 3.58 4.44 3.44 3.01 2.86 2.44];
CD0LT = [0.013 0.013 0.013 0.014 0.031 0.041 0.039 0.036 0.035];
etaLT = [0.54 0.54 0.54 0.75 0.79 0.78 0.89 0.93 0.93];

% Lookup Table 3: Propulsion Data (Thrust as function of mach number and altitude)
MachLT = [0; 0.2; 0.4; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8];
AltLT = 304.8*[0 5 10 15 20 25 30 40 50 70];
ThrustLT = 4448.2*[24.2 24.0 20.3 17.3 14.5 12.2 10.2 5.7 3.4 0.1;
28.0 24.6 21.1 18.1 15.2 12.8 10.7 6.5 3.9 0.2;
28.3 25.2 21.9 18.7 15.9 13.4 11.2 7.3 4.4 0.4;
30.8 27.2 23.8 20.5 17.3 14.7 12.3 8.1 4.9 0.8;
34.5 30.3 26.6 23.2 19.8 16.8 14.1 9.4 5.6 1.1;
37.9 34.3 30.4 26.8 23.3 19.8 16.8 11.2 6.8 1.4;
36.1 38.0 34.9 31.3 27.3 23.6 20.1 13.4 8.3 1.7;
36.1 36.6 38.5 36.1 31.6 28.1 24.2 16.2 10.0 2.2;
36.1 35.2 42.1 38.7 35.7 32.0 28.1 19.3 11.9 2.9;
36.1 33.8 45.7 41.3 39.8 34.6 31.1 21.7 13.3 3.1];

% Fitting of Aerodynamic data with piecewise splines
Clalphadat = pchip(MachLTAero,ClalphaLT);
CDdat = pchip(MachLTAero,CD0LT);
etadat = pchip(MachLTAero,etaLT);

Atomsrho=pchip(AtomsData(:,1),AtomsData(:,2));
Atomssos=pchip(AtomsData(:,1),AtomsData(:,3));


% Boundary Conditions 
alt0 = 0; 
altf = 19994.88; 
speed0 = 129.314; 
speedf = 295.092; 
fpa0 = 0; 
fpaf = 0; 
mass0 = 19050.864;

% Simple Bounds
altmin = 0; altmax = 21031.2;
speedmin = 5; speedmax = 1000;
fpamin = -40*pi/180; fpamax = 40*pi/180;
massmin = 22; massmax = 20410;
alphamin = -pi/4; alphamax = pi/4;

%%


%Initial Time. t0<tf
problem.time.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=0;     
problem.time.tf_max=400; 
guess.tf=324;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[alt0 speed0 fpa0 mass0];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[alt0 speed0 fpa0 mass0]; 
problem.states.x0u=[alt0 speed0 fpa0 mass0]; 

% State bounds. xl=< x <=xu
problem.states.xl=[altmin speedmin fpamin massmin]; 
problem.states.xu=[altmax speedmax fpamax massmax]; 

% State error bounds
problem.states.xErrorTol=[0.1 0.1 deg2rad(0.1) 0.1];

% State constraint error bounds
problem.states.xConstraintTol=[0.1 0.1 deg2rad(0.1) 0.1];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[altf speedf fpaf massmin]; 
problem.states.xfu=[altf speedf fpaf massmax];

% Guess the state trajectories with [x0 xf]
guess.time=[];
guess.states(:,1)=[alt0 altf];
guess.states(:,2)=[speed0 speedf];
guess.states(:,3)=[fpa0 fpaf];
guess.states(:,4)=[mass0 mass0];

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

% Input constraint error bounds
problem.inputs.uConstraintTol=[deg2rad(0.1)];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[20 -20]*pi/180;


% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[];
problem.constraints.gu=[];
problem.constraints.gTol=[];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];

% store the necessary problem parameters used in the functions
problem.data.CDdat = CDdat;
problem.data.Clalphadat = Clalphadat;
problem.data.etadat = etadat;
problem.data.M = MachLT;
problem.data.M2 = MachLTAero;
problem.data.alt = AltLT;
problem.data.T = ThrustLT;
problem.data.Re = 6378145;
problem.data.mu = 3.986e14;
problem.data.S = 49.2386;
problem.data.g0 = 9.80665;
problem.data.Isp = 1600;
problem.data.H = 7254.24;
problem.data.rho0 = 1.225;
problem.data.Atomsrho = Atomsrho;
problem.data.Atomssos = Atomssos;
[aa,mm] = meshgrid(AltLT,MachLT);
problem.data.aa = aa;
problem.data.mm = mm;

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc_unscaled,@b_unscaled};
problem.constraintErrorTol=[problem.constraints.gTol,problem.constraints.gTol,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,data)


% L_unscaled - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost=L_unscaled(x,xr,u,ur,p,t,data)
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

boundaryCost=-xf(4);

%------------- END OF CODE --------------


function dx = f_unscaled(x,u,p,t,data)

% f_unscaled - Returns the ODE right hand side where x'= f(x,u,p,t)
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = f_unscaled(x,u,p,t,data)
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

Atomsrho = data.Atomsrho;
Atomssos = data.Atomssos;
TLT = data.T;
mu = data.mu;
S = data.S;
g0 = data.g0;
Isp = data.Isp;
Re = data.Re;

h = x(:,1);
v = x(:,2);
fpa = x(:,3);
mass = x(:,4);
alpha = u(:,1);

r = h+Re;
rho = ppval(Atomsrho,h);
sos = ppval(Atomssos,h);
Mach = v./sos;

ii = Mach>=0.8;
jj = Mach<0.8;
mpoly = Mach(ii);
CD0 = zeros(length(Mach),1);
Clalpha = zeros(length(Mach),1);
eta = zeros(length(Mach),1);
if any(ii)
CD0(ii) = ppval(data.CDdat,mpoly);
Clalpha(ii) = ppval(data.Clalphadat,mpoly);
eta(ii) = ppval(data.etadat,mpoly);
end
if any(jj)
CD0(jj) = 0.013;
Clalpha(jj) = 3.44;
eta(jj) = 0.54;
end

Thrust = interp2(data.aa,data.mm,TLT,h,Mach,'spline');
CD = CD0 + eta.*Clalpha.*alpha.^2;
CL = Clalpha.*alpha;
q = 0.5.*rho.*v.*v;
D = q.*S.*CD;
L = q.*S.*CL;

hdot = v.*sin(fpa);
vdot = (Thrust.*cos(alpha)-D)./mass - mu.*sin(fpa)./r.^2;
fpadot = (Thrust.*sin(alpha)+L)./(mass.*v)+cos(fpa).*(v./r-mu./(v.*r.^2));
mdot = -Thrust./(g0.*Isp);

dx = [hdot, vdot, fpadot, mdot];

%------------- END OF CODE --------------


function c=g_unscaled(x,u,p,t,data)

% g_unscaled - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  c=g_unscaled(x,u,p,t,data)
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

c=[];

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

function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
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

if isfield(vdat,'gFilter')
    c(:,vdat.gFilter)=[];
end

if strcmp(vdat.mode.currentMode,'Feasibility')
    c=[c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np) c+p(:,end-vdat.mode.np+1:end)];
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