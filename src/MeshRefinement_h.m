function [ options ] = MeshRefinement_h( options, data, solution)
% MeshRefinement_h - Conduct mesh refinement for h typed discretization
% methods
%
%
% Syntax:  [ options ] = MeshRefinement_h( options, data, solution,minItervalScale )
%
% Inputs:
%    options - user supplied settings
%    data - Structure of data required to compute
%    solution - Data structure containing the solution
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


Error=solution.ErrorScaled;

% Intervals from the old mesh
if (strcmp(options.transcription,'hermite'))
    tau_inc=[0;data.tau_inc(2:2:end)];
else
    tau_inc=[0;data.tau_inc];
end

% Determine the minimum interval allowed for the new mesh
if min(diff(tau_inc))>1e-8
    mintimeinterval=min(diff(tau_inc))/2;
else
    mintimeinterval=min(diff(tau_inc));
end
mintimeinterval=max(options.mintimeinterval/solution.tf,mintimeinterval);

% Identify location of mesh refinement based on discretization errors
Mold=length(options.tau);
n=data.sizes{3};
idx=[];

n_addIntv=length(tau_inc); % maximum number of mesh point additions
for i=1:n
    [~,idxi]=sort(Error(:,i),'descend');
    idx_satis=Error(idxi,i)<options.discErrorTol_Scaled(i);
    idxi(idx_satis)=[];

    if length(idxi)>n_addIntv
        idxi=idxi(1:n_addIntv);
    end
    
    idx=union(idx,idxi);
end

% Also refine the neighbouring intervals if possible
if length(idx)<n_addIntv/2
    idx=unique([idx,idx-1,idx-2,idx+1,idx+2]);
    idx=idx(idx>0);
    idx=idx(idx<=size(Error,1));
end

% Identify refinement locations based on constraint violation
ConstraintError=solution.ConstraintError;
idxc=[];
for i=1:length(options.constraintErrorTol)
    [~,idxci] = sort(ConstraintError(:,i),'descend');
    idx_satis=ConstraintError(idxci,i)<=options.constraintErrorTol(i);
    idxci(idx_satis)=[];
    idxc=union(idxc,idxci);
end
if (strcmp(options.transcription,'hermite'))
    tau_remesh_constraint=solution.T_ConstraintError(idxc) / solution.tf*2;
else
    tau_remesh_constraint=solution.T_ConstraintError(idxc) / solution.tf;
end


% Refinement (subdividing intervals) based on discretization error
k=1;i=1;
while k<=n_addIntv && i<=length(idx)
    tau_remesh=data.tau_inc(idx(i));
    tau_mininterval= min(diff(data.tau_inc(max(1,idx(i)-1):min(length(Error(:,1)),idx(i)+1))))/2;
    if tau_mininterval<1e-8
        tau_mininterval=1e-8;
    end
    seg_remesh=tau_inc<tau_remesh;
    tau_front=tau_inc(seg_remesh);
    tau_rear=tau_inc(~seg_remesh);
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
        tau_inc=[tau_front;tau_remesh;tau_rear];
        k=k+1;i=1+1;
    else
        i=i+1;
    end
end

% Refinement (adding mesh point at exact location) based on constraint violation error
for j=1:length(tau_remesh_constraint)
    tau_remesh= tau_remesh_constraint(j);
    tau_mininterval= mintimeinterval;
    seg_remesh=tau_inc<tau_remesh;
    tau_front=tau_inc(seg_remesh);
    tau_rear=tau_inc(~seg_remesh);
    if ~isempty(tau_front) && tau_remesh-tau_front(end)<tau_mininterval
        if ~isempty(tau_rear) && tau_rear(1)-(tau_front(end)+tau_mininterval)>=tau_mininterval
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
        tau_inc=[tau_front;tau_remesh;tau_rear];
    end
end

% New mesh in the format of options structures
if (strcmp(options.transcription,'hermite'))
   options.tau=diff(tau_inc)/2;
else
   options.tau=diff(tau_inc);
end
        
% Increase the maximum number of iterations allowed accordingly
options.ipopt.max_iter=options.ipopt.max_iter*(length(options.tau))/(Mold);
end

