function [ const ] = userFunction_Adigator_LGR_All_Scaling( X_in,U,p,t,data)
%userFunction_Adigator_LGR_All_Scaling - Adigator template for user defined function (dynamics + eqaulity + inquality path constraints) with direct collocation method (hp-type) and automatic scaling
%
% Syntax:  [ const ] = userFunction_Adigator_LGR_All_Scaling( X_in,U,p,t,data)
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

f=data.data.InternalDynamics;



vdat=data.data;

X_Np1=X_in;
X=X_Np1(1:end-1,:);

dimX1=size(X,1);dimX2=size(X,2);
dimU1=size(U,1);dimU2=size(U,2);
XshiftMat=data.scaling.XshiftMat;
UshiftMat=data.scaling.UshiftMat;
X=data.scaling.XunscaleMat*X(:)-XshiftMat(:);
U=data.scaling.UunscaleMat*U(:)-UshiftMat(:);
% X=data.scaling.XunscaleMat*(X(:)-XshiftMat(:));
% U=data.scaling.UunscaleMat*(U(:)-UshiftMat(:));
X=reshape(X,dimX1,dimX2);
U=reshape(U,dimU1,dimU2);
if isfield(data.data,'Pscale')
%     p=(p-data.data.Pshift)./data.data.Pscale;
    p=p./data.data.Pscale-data.data.Pshift;
end

n=data.sizes{3};
ng=data.sizes{5};
M=data.sizes{7};
D_structure=data.map.D_structure;

t_0=t(1);
t_f=t(end);
delta_t=t_f-t_0;
T=delta_t/2*data.tau_inc+delta_t/2;
P=repmat(p,M,1);

[dynF_org,eqConst,neqConst]=f(X,U,P,T,vdat);
dynF=data.scaling.XscaleMat*dynF_org(:);
dynF=reshape(dynF,dimX1,dimX2);

g_all=reshape([eqConst;neqConst],M*ng,1);
g_vect=g_all(data.gAllidx);
t_segment=(t_f-t_0)/2*data.tau_segment'+(t_f+t_0)/2; %Time at start/end of each segment
t_segment_end=data.t_segment_mat_m*t_segment;
diag_t_segment_end=sparse(diag(t_segment_end));
const=[reshape(D_structure*X_Np1-diag_t_segment_end*dynF,M*n,1);
g_vect];


end

