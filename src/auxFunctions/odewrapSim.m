function dx=odewrapSim(t,x,solution,problem)
%odewrapSim - warp function for simulation of the obtained solution with ODE solvers
%
% Syntax:  dx=odewrapSim(t,x,solution,problem)
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


u=zeros(length(t),length(problem.sim.inputX)+length(problem.sim.inputU));
for i=1:length(problem.sim.inputX)
    u(:,i)=speval(solution,'X',i,t);
end
for i=1:length(problem.sim.inputU)
    u(:,i+length(problem.sim.inputX))=speval(solution,'U',i,t);
end

ul=[problem.states.xl(problem.sim.inputX) problem.inputs.ul(problem.sim.inputU)];
uu=[problem.states.xu(problem.sim.inputX) problem.inputs.uu(problem.sim.inputU)];
if any(u<ul)
    u(u<ul)=ul(u<ul);
end
if any(u>uu)
    u(u>uu)=uu(u>uu);
end

% Evaluate ODE right-hand side
p=repmat(solution.p,length(t),1);

f=problem.sim.functions;
dx=f(x',u,p,t,problem.data)';

end