function [ data, tau ] = genTimeMesh( options, data, ns, M )
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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


if options.tau==0
    tau=ns*ones(M-1,1)/(M-1); 
else
   texst=ones(ns,1);
   tau=kron(options.tau,texst);
end
data.tau_inc=cumsum(tau);
if abs(sum(tau)-ns)>sqrt(eps);error('Time vector (tau) should sum to 1');end
data.tau=tau;

end

