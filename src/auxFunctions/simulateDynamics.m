function [ x_crt ] = simulateDynamics( problem, simtime, u_crt, x_crt, tstep, odesolver )
% simulateSolution - Simulate the obtained solution (open-loop) with ODE integration
%
% Syntax:  
%          [ tv, xv, uv ] = simulateSolution( problem, solution)	
%          [ tv, xv, uv ] = simulateSolution( problem, solution, odesolver)                     (Request a specific ODE solver)
%          [ tv, xv, uv ] = simulateSolution( problem, solution, odesolver, tstep)              (Request a specific ODE solver and integration time step)
%          [ tv, xv, uv ] = simulateSolution( problem, solution, odesolver, tstep, s_i_swarp)   (Request a specific ODE solver and integration time step, with state-input swarping)
% 
% Output:
%    tv - time vector 
%    xv - state vector
%    uv - input vector
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
T=[ simtime simtime+tstep];

ul=problem.inputs.ul(problem.sim.inputU);
uu=problem.inputs.uu(problem.sim.inputU);
if u_crt<ul
    u_crt=ul;
elseif u_crt>uu
    u_crt=uu;
end
    
if strcmp(odesolver,'ode113')
    [~,xv]=ode113(@(t,x)odewrapSimDynamics(t,x,u_crt,problem),T,x_crt);
elseif strcmp(odesolver,'ode45') 
    [~,xv]=ode45(@(t,x)odewrapSimDynamics(t,x,u_crt,problem),T,x_crt);
elseif strcmp(odesolver,'ode23') 
    [~,xv]=ode23(@(t,x)odewrapSimDynamics(t,x,u_crt,problem),T,x_crt);
else
    error('Select ODE solver not supported');
end

x_crt=xv(end,:);

end

