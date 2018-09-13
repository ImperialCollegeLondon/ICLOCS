function [ data, tau_inc, tau_seg, tau, LGR ] = genTimeMeshLGR( options, data, nps, npdu, npd, npduidx, M )
%genTimeMesh - generate discretized mesh along the time dimension, for p/hp methods
%
% Syntax:  [ data, tau_inc, tau_seg, tau, LGR ] = genTimeMeshLGR( options, data, nps, npdu, npd, npduidx, M )
%
% Inputs:
%    options, data, nps, npdu, npd, npduidx, M - Defined in transcribeOCP
%
% Outputs:
%    data - Defined in transcribeOCP
%    tau_inc - time vector representing the time mesh, in cumulative sum format 
%    tau_seg - interval end points
%    tau - time vector representing interval lengths of mesh, should sum to 1
%    LGR - a structure containing the information fo the LGR points,
%    weights, differention matrix, etc
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


if isfield(options,'tau_segment')
    tau_segment=options.tau_segment; %user specified mesh segments
else
    tau_segment=-1:2/nps:1; %uniform mesh segments
end
%Calculate all the LGR points, weights and differentiation matrices
LGR=generateLGR_All(npdu); 
tau_local_seg=zeros(M,1); %mesh in local (segmental) time frame
idxtemp=[0,cumsum(npd)];
for i=1:length(npduidx)
    tau_local_seg(idxtemp(i)+1:idxtemp(i+1),:)=LGR.points{npduidx(i)};
end
tau_inc=zeros(M,1); %mesh in incremental (global) time frame
tau_seg=cell(length(tau_segment)-1,1);
for i=1:length(tau_segment)-1
    tau_seg{i} = ((tau_segment(i+1)+tau_segment(i)) + (tau_segment(i+1)-tau_segment(i)) * [LGR.points{npduidx(i)}; 1]) / (1-(-1));
    tau_inc(idxtemp(i)+1:idxtemp(i+1),:)=tau_seg{i}(1:end-1);
end
tau=diff(tau_inc); %mesh in difference formulation
data.tau_inc=tau_inc;
data.tau_local_seg=tau_local_seg;
data.tau_segment=tau_segment;

data.tau=tau;

end

