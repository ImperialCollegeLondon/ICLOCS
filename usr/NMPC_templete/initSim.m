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
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

clear all

addpath('MPCController1')
addpath('MPCController2')
...

global infoNLP data simtime solution tstep N_node
simtime=...; % initialize simulation time
tstep=...; % time step of simulation
N_node=...; % initial number of mesh points

options= settings_myProblem_Phase1(N_node);                  % Get options and solver settings
[problem,guess]=myProblem_Phase1;          % Fetch the problem definition

% initialize problem
solution.tf=problem.time.tf_max; 
solution.status.status=0;
solution.t_ref=0;
par=problem.data;
[infoNLP,data]=transcribeOCP(problem,guess,options); 
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});