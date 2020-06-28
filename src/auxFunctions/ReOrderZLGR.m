function [ data ] = ReOrderZLGR( data )
%ReOrderZLGR - Generate alternative ordering for LGR direct transcription
%formulation
%
% Syntax:  data = ReOrderZLGR( data )
%
% Inputs:
%    data - The data structure before reordering
%
% Outputs:
%    data - The data structure after reordering
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

[nt,np,n,m,ng,nb,M,~,~,~,~,~,~,nrcl,nrcu,nrce,~]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;

z_idx_org=transpose(1:((M+1)*n+M*m+np+nt));

X_Np1=reshape(data.map.Vx*z_idx_org,M+1,n);
U=reshape(data.map.Vu*z_idx_org,M,m);

vert_idx_org=transpose(1:(M*(n+ng)+nb+nrc));
vert_idx_rsp=reshape(vert_idx_org(1:M*(n+ng)),M,(n+ng));

z_idx_new=zeros(1,(m+n)*M);
vert_idx_new=zeros(1,M*(n+ng));
for i=1:M
    z_idx_new(1,(i-1)*(m+n)+1:i*(m+n))=[X_Np1(i,:),U(i,:)];
    vert_idx_new(1,(i-1)*(n+ng)+1:i*(n+ng))=vert_idx_rsp(i,:);
end

z_idx_new=[z_idx_new,X_Np1(end,:),z_idx_org(end-nt-np+1:end)']';
vert_idx_new=[vert_idx_new,vert_idx_org(M*(n+ng)+1:end)']';
z_idx_back=arrayfun(@(x)find(z_idx_new==x,1),z_idx_org);
vert_idx_back=arrayfun(@(x)find(vert_idx_new==x,1),vert_idx_org);


data.reorder.z_idx=z_idx_new;
data.reorder.z_idx_back=z_idx_back;
data.reorder.vert_idx=vert_idx_new;
data.reorder.vert_idx_back=vert_idx_back;
end

