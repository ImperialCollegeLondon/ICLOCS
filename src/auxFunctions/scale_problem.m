function [ problem,guess ] = scale_problem( problem,guess )
%scale_problem - automatic scaling of the problem
%
% Syntax:  [ problem,guess ] = scale_problem( problem,guess )
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

if isfield(problem,'scaling')
     xl=problem.scaling.states.lb;
     xu=problem.scaling.states.ub;

     idxeq=xl==xu;
     xl(idxeq)=xl(idxeq)-0.1;
     xu(idxeq)=xu(idxeq)+0.1;

     problem.states.scales=1./(xu-xl);
     problem.states.shifts=0.5-xu./(xu-xl);
    
elseif ~isfield(problem.states,'scales')
    
     xl=problem.states.xl;xl(isinf(xl))=-500;
     xu=problem.states.xu;xu(isinf(xu))=500;

     idxeq=xl==xu;
     xl(idxeq)=xl(idxeq)-0.1;
     xu(idxeq)=xu(idxeq)+0.1;

     problem.states.scales=1./(xu-xl);
     problem.states.shifts=0.5-xu./(xu-xl);
end
 
 if isfield(problem.states,'x0') && ~isempty(problem.states.x0)
    problem.states.x0 = scale_variables( problem.states.x0, problem.states.scales, problem.states.shifts );
 end
 problem.states.x0l = scale_variables( problem.states.x0l, problem.states.scales, problem.states.shifts );
 problem.states.x0u = scale_variables( problem.states.x0u, problem.states.scales, problem.states.shifts );
 problem.states.xl = scale_variables( problem.states.xl, problem.states.scales, problem.states.shifts );
 problem.states.xu = scale_variables( problem.states.xu, problem.states.scales, problem.states.shifts );
 problem.states.xfl = scale_variables( problem.states.xfl, problem.states.scales, problem.states.shifts );
 problem.states.xfu = scale_variables( problem.states.xfu, problem.states.scales, problem.states.shifts );
 if ~isempty(guess.states)
    guess.states=scale_variables(guess.states, problem.states.scales, problem.states.shifts );
 end
 if isfield(problem.states,'xrl')
     problem.states.xrl=scale_variables( problem.states.xrl, problem.states.scales, 0 );
     problem.states.xru=scale_variables( problem.states.xru, problem.states.scales, 0 );
 end
 
 if ~isempty(problem.setpoints.states)
    problem.setpoints.states=scale_variables( problem.setpoints.states, problem.states.scales, problem.states.shifts );
 end

 if isfield(problem,'scaling')
     uu=problem.scaling.inputs.ub;
     ul=problem.scaling.inputs.lb;

     idxeq=ul==uu;
     ul(idxeq)=ul(idxeq)-0.1;
     uu(idxeq)=uu(idxeq)+0.1;

     problem.inputs.scales=1./(uu-ul);
     problem.inputs.shifts=0.5-uu./(uu-ul);
     
 elseif ~isfield(problem.inputs,'scales')
     uu=problem.inputs.uu;uu(isinf(uu))=500;
     ul=problem.inputs.ul;ul(isinf(ul))=-500;

     idxeq=ul==uu;
     ul(idxeq)=ul(idxeq)-0.1;
     uu(idxeq)=uu(idxeq)+0.1;

     problem.inputs.scales=1./(uu-ul);
     problem.inputs.shifts=0.5-uu./(uu-ul);
 end

 problem.inputs.ul=scale_variables( problem.inputs.ul, problem.inputs.scales, problem.inputs.shifts );
 problem.inputs.uu=scale_variables( problem.inputs.uu, problem.inputs.scales, problem.inputs.shifts );
 problem.inputs.u0l=scale_variables( problem.inputs.u0l, problem.inputs.scales, problem.inputs.shifts );
 problem.inputs.u0u=scale_variables( problem.inputs.u0u, problem.inputs.scales, problem.inputs.shifts );
 
 if isfield(problem.inputs,'url')
     problem.inputs.url=scale_variables( problem.inputs.url, problem.inputs.scales, 0 );
     problem.inputs.uru=scale_variables( problem.inputs.uru, problem.inputs.scales, 0 );
 end
 if ~isempty(guess.inputs)
    guess.inputs=scale_variables( guess.inputs, problem.inputs.scales, problem.inputs.shifts );
 end
 if ~isempty(problem.setpoints.inputs)
    problem.setpoints.inputs=scale_variables(problem.setpoints.inputs, problem.inputs.scales, problem.inputs.shifts );
 end

 if ~isempty(guess.parameters)
     
     if ~isfield(problem.parameters,'scales')
         pu=problem.parameters.pu;pu(isinf(pu))=500;
         pl=problem.parameters.pl;pl(isinf(pl))=-500;

         idxeq=pl==pu;
         pl(idxeq)=pl(idxeq)-0.1;
         pu(idxeq)=pu(idxeq)+0.1;

         problem.parameters.scales=1./(pu-pl);
         problem.parameters.shifts=0.5-pu./(pu-pl);
     end
     problem.parameters.pu = scale_variables( problem.parameters.pu, problem.parameters.scales, problem.parameters.shifts );
     problem.parameters.pl = scale_variables( problem.parameters.pl, problem.parameters.scales, problem.parameters.shifts );
     guess.parameters = scale_variables( guess.parameters, problem.parameters.scales, problem.parameters.shifts );
 end
     
     
 
%  problem.time.scales=1./problem.time.tf_max;
%  problem.time.scales=1./problem.time.tf_max*200;
%  problem.time.shifts=0.5-problem.time.tf_max./(problem.time.tf_max-problem.time.tf_min);
 problem.time.scales=1;
 problem.time.shifts=0;

 problem.time.tf_max=scale_variables( problem.time.tf_max, problem.time.scales, problem.time.shifts );
 problem.time.tf_min=scale_variables( problem.time.tf_min, problem.time.scales, problem.time.shifts );
 guess.tf=scale_variables( guess.tf, problem.time.scales, problem.time.shifts );
%  
end

