function [ tv, xv, uv ] = simulateSolution( varargin )
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

if nargin==2
    problem=varargin{1};
    solution=varargin{2};
    T=[solution.t0 solution.tf];
    x0=solution.x0;
elseif nargin==3
    problem=varargin{1};
    solution=varargin{2};
    odesolver=varargin{3};
    T=[solution.t0 solution.tf];
    x0=solution.x0;
elseif nargin==4
    problem=varargin{1};
    solution=varargin{2};
    odesolver=varargin{3};
    tstep=varargin{4};
    T=solution.t0:tstep:solution.tf;
    x0=solution.x0;
elseif nargin==5
    problem=varargin{1};
    solution=varargin{2};
    odesolver=varargin{3};
    tstep=varargin{4};
    T=solution.t0:tstep:solution.tf;
    s_i_swarp=varargin{5};
    state_rem=s_i_swarp(1,:);
    input_rem=s_i_swarp(2,:);
    solution.Up(input_rem)=solution.Xp(state_rem);
    solution.Xp(state_rem)=[];
    x0=solution.x0;
    x0(state_rem)=[];
else
    error('Number of input parameter to the simulation function not vaild');
end
    
if nargin<3 || strcmp(odesolver,'ode113') 
    [tv,xv]=ode113(@(t,x)odewrapSim(t,x,solution,problem),T,x0);
elseif strcmp(odesolver,'ode45') 
    [tv,xv]=ode45(@(t,x)odewrapSim(t,x,solution,problem),T,x0);
elseif strcmp(odesolver,'ode23') 
    [tv,xv]=ode23(@(t,x)odewrapSim(t,x,solution,problem),T,x0);
else
    error('Select ODE solver not supported');
end

uv=zeros(length(tv),length(problem.sim.inputX)+length(problem.sim.inputU));
for i=1:length(problem.sim.inputX)
    uv(:,i)=speval(solution,'X',problem.sim.inputX(i),tv);
end
for i=1:length(problem.sim.inputU)
    uv(:,i+length(problem.sim.inputX))=speval(solution,'U',problem.sim.inputU(i),tv);
end

ul=[problem.states.xl(problem.sim.inputX) problem.inputs.ul(problem.sim.inputU)];
uu=[problem.states.xu(problem.sim.inputX) problem.inputs.uu(problem.sim.inputU)];
if any(uv<ul)
    ul_mat=repmat(ul,size(uv,1));
    uv(uv<ul)=ul_mat(uv<ul);
end
if any(uv>uu)
    uu_mat=repmat(uu,size(uv,1));
    uv(uv>uu)=uu_mat(uv>uu);
end

end

