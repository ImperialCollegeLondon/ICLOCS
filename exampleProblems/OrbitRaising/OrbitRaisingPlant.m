function dx = OrbitRaisingPlant(x,u,p,t,data)

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

r=x(:,1);theta=x(:,2);v_r=x(:,3);v_theta=x(:,4);
u_1=u(:,1);u_2=u(:,2);

T1=data.T1;
md=data.md;
dx=[v_r,v_theta./r,(v_theta.*v_theta)./r-r.^(-2)+(T1./(1-md*t)).*u_1,...
    (-v_r.*v_theta)./r+(T1./(1-md*t)).*u_2];

