function [sys,x0,str,ts] = CartPoleSwingUp_sfun(t,x,u,flag,varargin)


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
   [sys,DERX] = CartPoleSwingUp_eqn(t,x,u,varargin);
   
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
sizes.NumInputs      = 1;
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

function [Y,DERX] = CartPoleSwingUp_eqn(t,x,u,varargin);

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
    vdat = varargin{1}{1};    
end

x1 = x(1);x2 = x(2);x3 = x(3);x4 = x(4);
u1 = u(1);

% Model mismatch
% L=vdat.L;
% m1 = vdat.m1;
% m2 = vdat.m2;
% g=vdat.g;

L=vdat.L+0.05;
m1 = vdat.m1+0.1;
m2 = vdat.m2-0.02;
g=vdat.g;

dx1 = (L.*m2.*sin(x4).*x2.^2 + u1 + m2.*g.*cos(x4).*sin(x4))./...
            (m1+m2.*(1-cos(x4).^2));

dx2 = -( L.*m2.*cos(x4).*sin(x4).*x2.^2 + u1.*cos(x4) + (m1+m2).*g.*sin(x4))./...
            (L.*m1+L.*m2.*(1-cos(x4).^2));

dx3 = x1;

dx4 = x2;

DERX = [dx1, dx2, dx3, dx4];
Y = [x1,x2,x3,x4];