function [ options ] = MeshRefinement_Auto( options, data, solution, minItervalScale )
% MeshRefinement_Auto - Conduct mesh refinement for h typed discretization
% methods
%
%
% Syntax:  [ options ] = MeshRefinement_h( options, data, solution,minItervalScale )
%
% Inputs:
%    options - user supplied settings
%    data - Structure of data required to compute
%    solution - Data structure containing the solution
%    minItervalScale - scaling factor for the minimum mesh size allowed
%    
% Outputs:
%    options - settings with the new (refined) mesh
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

Error=solution.Error;

options.tau_segment=data.tau_segment;
options.npsegment=data.sizes{1, 10};

n_addIntv=6;
npd_Max=15;
npd_Min=3;

% Intervals from the old mesh
tau_seg_idx=find(data.tau_local_seg==-1);

n=data.sizes{3};
idx=[];

% Identify location of mesh refinement based on discretization errors
for i=1:n
    [~,idxi]=sort(Error(:,i),'descend');
    idx_satis=Error(idxi,i)<options.discErrorTol_Scaled(i);
    idxi(idx_satis)=[];

    if length(idxi)>n_addIntv
        idxi=idxi(1:n_addIntv);
    end
    
    idx=union(idx,idxi);
end
T_remesh_div=solution.T(idx+1);

% Identify refinement locations based on constraint violation
ConstraintError=solution.ConstraintError;
idxc=[];
for i=1:length(options.constraintErrorTol)
    [~,idxci] = findpeaks(ConstraintError(:,i),'SortStr','descend');
    idx_satis=ConstraintError(idxci,i)<options.constraintErrorTol(i);
    idxci(idx_satis)=[];
    idxc=union(idxc,idxci);
end
T_remesh_constraint=solution.T_ConstraintError(idxc);

% Identify location of mesh refinement
mesh_AI_indicator=zeros(1,length(options.tau_segment)-1);
for i=1:length(idx)
    for j=1:length(tau_seg_idx)
        if idx(i)>=tau_seg_idx(j) && j~=length(tau_seg_idx) && idx(i)<tau_seg_idx(j+1) 
            mesh_AI_indicator(j)=1;
        elseif idx(i)>=tau_seg_idx(j) && j==length(tau_seg_idx)
            mesh_AI_indicator(j)=1;
        end
    end      
end

Max_error=zeros(length(options.tau_segment)-1,n);
for i=1:length(options.tau_segment)-1
    if i~=length(options.tau_segment)-1
        Max_error(i,:)=max(Error(tau_seg_idx(i):(tau_seg_idx(i+1)-1),:),[],1);
    else
        Max_error(i,:)=max(Error(tau_seg_idx(i):end,:),[],1);
    end
    % increase polynomial order for intervals without mesh refinement, if
    % discretization error still large
    if any(Max_error(i,:)>options.discErrorTol) && mesh_AI_indicator(i)~=1
        T_remesh=solution.T(1+tau_seg_idx(i));
        tau_remesh= ((1-(-1))*T_remesh-(solution.T(1)+solution.tf))/(solution.tf-solution.T(1));
        seg_remesh=options.tau_segment<tau_remesh;
        segidx_remesh=length(options.tau_segment(seg_remesh));
        p_inc=ceil(max(log10(Max_error(i,:)./options.discErrorTol)));
            if options.npsegment(segidx_remesh)+p_inc<=npd_Max
                options.npsegment(segidx_remesh)=options.npsegment(segidx_remesh)+p_inc;
            elseif options.npsegment(segidx_remesh)<=npd_Max
                options.npsegment(segidx_remesh)=npd_Max;
            else
                T_remesh_div=union(T_remesh_div,(solution.T(tau_seg_idx(i))+solution.T(1+tau_seg_idx(i)))/2);
            end

            if isfield(options, 'pdegree')
                options = rmfield(options,{'pdegree','nsegment'});
            end

            if options.adaptseg
                options.nsegment=length(options.npsegment);
            end
    % decrease polynomial order for intervals without mesh refinement, if
    % discretization error smaller than required
    elseif all(Max_error(i,:)<options.discErrorTol) && mesh_AI_indicator(i)~=1
        T_remesh=solution.T(1+tau_seg_idx(i));
        tau_remesh= ((1-(-1))*T_remesh-(solution.T(1)+solution.tf))/(solution.tf-solution.T(1));
        seg_remesh=options.tau_segment<tau_remesh;
        segidx_remesh=length(options.tau_segment(seg_remesh));
        p_dec=ceil(max(log10(sqrt(Max_error(i,:)./options.discErrorTol))));
            if options.npsegment(segidx_remesh)+p_dec>=npd_Min
                options.npsegment(segidx_remesh)=options.npsegment(segidx_remesh)+p_dec;
            elseif options.npsegment(segidx_remesh)>npd_Min
                options.npsegment(segidx_remesh)=npd_Min;
            end

            if isfield(options, 'pdegree')
                options = rmfield(options,{'pdegree','nsegment'});
            end

            if options.adaptseg
                options.nsegment=length(options.npsegment);
            end
    end
end

% Refine the mesh where needed
T_remesh_div=[T_remesh_div;T_remesh_constraint'];
for i=1:length(T_remesh_div)
    T_remesh=T_remesh_div(i);
    tau_remesh= ((1-(-1))*T_remesh-(solution.T(1)+solution.tf))/(solution.tf-solution.T(1));
    tau_mininterval= options.mintimeinterval*minItervalScale/(solution.tf-solution.T(1))*(1-(-1));
    tau_mininterval=max(tau_mininterval,1e-9);
    seg_remesh=options.tau_segment<tau_remesh;
    tau_front=options.tau_segment(seg_remesh);
    tau_rear=options.tau_segment(~seg_remesh);
    if ~isempty(tau_front) && tau_remesh-tau_front(end)<tau_mininterval
        if ~isempty(tau_rear) && ~isempty(tau_front) && tau_rear(1)-(tau_front(end)+tau_mininterval)>=tau_mininterval
            tau_remesh=tau_front(end)+tau_mininterval;
        else
            tau_remesh=[];
        end
    elseif ~isempty(tau_rear) && tau_rear(1)-tau_remesh<tau_mininterval
        if ~isempty(tau_rear) && ~isempty(tau_front) && (tau_rear(1)-tau_mininterval)-tau_front(end)>=tau_mininterval
            tau_remesh=tau_rear(1)-tau_mininterval;
        else
            tau_remesh=[];
        end
    end

    if ~isempty(tau_remesh)
        segidx_remesh=length(options.tau_segment(seg_remesh));
        options.tau_segment=[tau_front tau_remesh tau_rear];
        new_pd=max(ceil(options.npsegment(segidx_remesh))/2,8);
        options.npsegment=[options.npsegment(1:segidx_remesh-1) new_pd new_pd options.npsegment(segidx_remesh+1:end)];

        if isfield(options, 'pdegree')
            options = rmfield(options,{'pdegree','nsegment'});
        end

        if options.adaptseg
            options.nsegment=length(options.npsegment);
        end
    end
end

end

