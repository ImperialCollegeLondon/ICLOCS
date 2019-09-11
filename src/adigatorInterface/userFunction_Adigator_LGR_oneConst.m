function [ const ] = userFunction_Adigator_LGR_oneConst( X_in,U,p,t,data)
%userFunction_Adigator_LGR_oneConst - Adigator template for user defined function (dynamics + equality or inequality path constraints) with direct collocation method (hp-type)
%
% Syntax:  [ const ] = userFunction_Adigator_LGR_oneConst( X_in,U,p,t,data)
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

n=data.sizes{3};
ng=data.sizes{5};
M=data.sizes{7};
D_structure=data.map.D_structure;

t_0=t(1);
t_f=t(end);
delta_t=t_f-t_0;
T=delta_t/2*data.tau_inc+delta_t/2;
P=repmat(p,M,1);

[dynF,gConst]=f(X,U,P,T,vdat);
g_all=reshape(gConst,M*ng,1);
g_vect=g_all(data.gAllidx);
t_segment=(t_f-t_0)/2*data.tau_segment'+(t_f+t_0)/2; %Time at start/end of each segment
t_segment_end=data.t_segment_mat_m*t_segment;
diag_t_segment_end=sparse(diag(t_segment_end));
const=[reshape(D_structure*X_Np1-diag_t_segment_end*dynF,M*n,1);
g_vect];


end

