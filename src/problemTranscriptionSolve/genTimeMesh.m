function [ data, tau ] = genTimeMesh( problem, options, data, ns, M )
%genTimeMesh - generate discretized mesh along the time dimension, for h
%methods
%
% Syntax:  [ data, tau ] = genTimeMesh( options, data, ns, M )
%
% Inputs:
%    options, data, ns, M - Defined in transcribeOCP
%
% Outputs:
%    data - Defined in transcribeOCP
%    tau - time vector representing interval lengths of mesh, should sum to 1
%     
%
% Other m-files required: none
% Subfunctions: none
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

%------------- BEGIN CODE --------------


if options.tau==0
    tau=ns*ones(M-1,1)/(M-1); 
    data.tau_inc=linspace(0,ns,M)';
    data.tau_inc(1)=[];
else
   texst=ones(ns,1);
   tau=kron(options.tau,texst);
   data.tau_inc=cumsum(tau);
end

if abs(sum(tau)-ns)>sqrt(eps);error('Time vector (tau) should sum to 1');end
data.tau=tau;


   if problem.time.t0_min==problem.time.t0_min && problem.time.tf_min==problem.time.tf_min
       if isfield(problem.inputs,'t_zone') && problem.inputs.t_zone>0
           t0=problem.time.t0_min;
           tf=problem.time.tf_min;
           T=(tf-t0)*[0;data.tau_inc]*data.Nm/ns+t0;
           [hc,edges,bin]  = histcounts(T,t0:problem.inputs.t_zone:tf);
           temp=cumsum([1 hc]);
           data.t_zone_map=temp(bin);       
       end
       data.ProblemTypes.FixedTime=1;
   else
       if isfield(problem.inputs,'t_zone')
           error('Zoning of input is only availabe for fixed t0 and tf problems')
       end
       data.ProblemTypes.FixedTime=0;
   end


end

