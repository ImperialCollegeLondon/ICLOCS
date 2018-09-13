function [sys,x0,str,ts] = MinEnergyClimb_sfun(t,x,u,flag,varargin)


% Keep DERX persistent. In this way, rcam_eqm has to be 
% evaluated once per simulation step.

persistent DERX 

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(varargin);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
   % sys=mdlDerivatives(t,x,u);
   sys = DERX;
   
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=[];

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
   % [sys]=mdlOutputs(t,x,u);
   [sys,DERX] = MinEnergyClimb_eqn(t,x,u,varargin);
   
  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=[];

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(varargin)

sizes = simsizes;

sizes.NumContStates  = 4;
sizes.NumDiscStates  =  0;
sizes.NumOutputs     = 4;
sizes.NumInputs      = 4;
sizes.DirFeedthrough =  1;
sizes.NumSampleTimes =  1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%

if nargin > 0,
  x0  = varargin{1}{2};
else
  x0  = [0;0;0;0];
end;

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
%%%%% function sys=mdlDerivatives(t,x,u)
%
%sys =   DERX;
%
%% end mdlDerivatives

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

%%%%% function [sys]=mdlOutputs(t,x,u,varargin)

% end mdlOutputs

% Instead, use rcam_eqm, which only has one more output!

function [Y,DERX] = MinEnergyClimb_eqn(t,x,u,varargin);

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
if nargin > 3
    data = varargin{1}{1};    
end

Atomsrho = data.Atomsrho;
Atomssos = data.Atomssos;
TLT = data.T;
mu = data.mu;
S = data.S;
g0 = data.g0;
Isp = data.Isp;
Re = data.Re;

h = x(1);
v = x(2);
fpa = x(3);
mass = x(4);
u_g= u(2);
w_g= u(3);
q_g= u(4);
alpha=u(1);
v=sqrt((v*sin(fpa)+w_g)^2+(v*cos(fpa)+u_g)^2);

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

DERX = [hdot, vdot, fpadot, mdot];
Y = [h,v,fpa,mass];