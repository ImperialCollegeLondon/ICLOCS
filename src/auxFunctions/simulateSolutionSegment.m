function [ tv, xv, uv ] = simulateSolution( varargin )
%simulateSolutionSegment Simulate a solution over an explicit time vector.
%
% Syntax:  
%   [tv,xv,uv] = simulateSolutionSegment(problem,solution,odesolver,T)
%   [tv,xv,uv] = simulateSolutionSegment(problem,solution,odesolver,T,stateInputSwap)
% 
% Inputs:
%   problem        - Original ICLOCS problem structure.
%   solution       - Solution structure returned by solveMyProblem.
%   odesolver      - ODE solver name: 'ode113', 'ode45', or 'ode23'.
%   T              - Explicit simulation time vector.
%   stateInputSwap - Optional two-row index map for examples that swap selected
%                    states and inputs during simulation.
%
% Outputs:
%   tv - Time vector returned by the ODE solver.
%   xv - Simulated state trajectory.
%   uv - Interpolated input trajectory evaluated at tv.
%
% Notes:
%   This function is the segment-time-vector variant of simulateSolution. It
%   simulates the supplied solution and does not re-optimize.
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
    T=varargin{4};
    x0=solution.x0;
elseif nargin==5
    problem=varargin{1};
    solution=varargin{2};
    odesolver=varargin{3};
    T=varargin{4};
    s_i_swarp=varargin{5};
    state_rem=s_i_swarp(1,:);
    input_rem=s_i_swarp(2,:);
    solution.Up(input_rem)=solution.Xp(state_rem);
    solution.Xp(state_rem)=[];
    x0=solution.x0;
    x0(state_rem)=[];
elseif nargin==6
    problem=varargin{1};
    solution=varargin{2};
    odesolver=varargin{3};
    T=varargin{4};
    s_i_swarp=varargin{6};
    state_add=s_i_swarp(2,:);
    input_add=s_i_swarp(1,:);
    solution.Xp(state_add)=solution.Up(input_add);
    solution.Up(input_add)=solution.dUp(input_add);
    x0=[solution.x0 solution.U(1,input_add)];
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
