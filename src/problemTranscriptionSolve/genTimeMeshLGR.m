function [ data, tau_inc, tau_seg, tau, LGR ] = genTimeMeshLGR( problem, options, data, nps, npdu, npd, npduidx, M )
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
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


if isfield(options,'tau_segment')
    tau_segment=options.tau_segment; %user specified mesh segments
else
    tau_segment=-1:2/nps:1; %uniform mesh segments
end
tau_ratio=(tau_segment+1)/2;
%Calculate all the LGR points, weights and differentiation matrices
LGR=generateLGR_All(npdu); 
tau_local_seg=zeros(M,1); %mesh in local (segmental) time frame
t_segment_mat_m=zeros(M,length(npduidx)+1);
t_segment_mat_p=zeros(M,length(npduidx)+1);
idxtemp=[0,cumsum(npd)];
for i=1:length(npduidx)
    tau_local_seg(idxtemp(i)+1:idxtemp(i+1),:)=LGR.points{npduidx(i)};
    tau_order_idx(idxtemp(i)+1:idxtemp(i+1),:)=transpose(1:npd(i));
    tau_segment_ratio_diff(idxtemp(i)+1:idxtemp(i+1),:)=tau_ratio(i+1)-tau_ratio(i);
    t_segment_mat_m(idxtemp(i)+1:idxtemp(i+1),i)=-0.5;
    t_segment_mat_m(idxtemp(i)+1:idxtemp(i+1),i+1)=0.5;
    t_segment_mat_p(idxtemp(i)+1:idxtemp(i+1),i)=0.5;
    t_segment_mat_p(idxtemp(i)+1:idxtemp(i+1),i+1)=0.5;
    if mod(i,3)==1
        tau_order_idx_hes(idxtemp(i)+1:idxtemp(i+1),:)=transpose(1:npd(i));
        if i==length(npduidx)
            tau_order_idx_hes(end+1,:)=101;
        end
    elseif mod(i,3)==2
        tau_order_idx_hes(idxtemp(i)+1:idxtemp(i+1),:)=transpose(101:100+npd(i));
        if i==length(npduidx)
            tau_order_idx_hes(end+1,:)=201;
        end
    else
        tau_order_idx_hes(idxtemp(i)+1:idxtemp(i+1),:)=transpose(201:200+npd(i));
        if i==length(npduidx)
            tau_order_idx_hes(end+1,:)=1;
        end
    end
    tau_seg_idx(idxtemp(i)+1:idxtemp(i+1),:)=i;
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
data.tau_order_idx=tau_order_idx;
data.tau_seg_idx=tau_seg_idx;
data.tau_order_idx_hes=tau_order_idx_hes;
data.tau=tau;
data.tau_segment_ratio_diff=tau_segment_ratio_diff;
data.t_segment_mat_m=sparse(t_segment_mat_m);
data.t_segment_mat_p=sparse(t_segment_mat_p);

if isfield(problem.inputs,'t_zone') && problem.inputs.t_zone>0
   if problem.time.t0_min==problem.time.t0_min && problem.time.tf_min==problem.time.tf_min
       t0=problem.time.t0_min;
       tf=problem.time.tf_min;
       T=(tf-t0)/2*tau_inc+(tf+t0)/2;
       [hc,edges,bin]  = histcounts(T,t0:problem.inputs.t_zone:tf);
       temp=cumsum([1 hc]);
       data.t_zone_map=temp(bin);
   else
       error('Zoning of input is only availabe for fixed t0 and tf problems')
   end
end

end

