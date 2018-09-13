%initSim - Initialization for Closed-loop Simulation with ICLOCS2
%
% Input:
%    tpx0 - vector in format [current time; current state]
%
% Output:
%    tfpu - vector in format [terminal time; computation time; input to implement]
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

clear all

global infoNLP data simtime solution tstep opt_step N opt_t_min
simtime=0; % initialize simulation time
tstep=0.01; % time step of simulation
opt_step=0; % time interval for re-optimization, set to zero for update ASAP
opt_t_min=4; % stop re-optimizing if solution.tf<opt_t_min
N=20; % initial number of mesh points

options= settings_h(N);                  % Get options and solver settings
[problem,guess]=MinEnergyClimbBryson;          % Fetch the problem definition

% initialize problem
solution.tf=problem.time.tf_max; 
solution.status.status=0;
solution.t_ref=0;
par=problem.data;
[infoNLP,data]=transcribeOCP(problem,guess,options); 
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});


