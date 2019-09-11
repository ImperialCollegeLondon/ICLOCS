function grad=gradientFD(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,t0,tf,data)
%GRADIENTFD - Generate gradient of the cost using finite-differences
%
% Syntax:  grad=gradientFD(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    Defined in directCollocation.m
%
% Outputs:
%    solution - Data structure containing the solution
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

e=data.options.perturbation.J;                                 % Pertubation size

[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});

nx=M*n;                               % Number of unknown states
nu=N*m;                               % Number of unknown controls
nz=nt+np+nx+nu;                       % Length of the optimization variable 
vdat=data.data;                       % Variables used in the functions 
t0=data.t0;
k0=data.k0;
% Compute contribution of L to the gradient

Lz=zeros(1,nz);
idx=data.FD.index.Ly;
nfd=size(idx,2);                               

et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;

%This vector is used to adjust the size of the vector of the derivative
%when the stage cost is identically zero 

snm=ones(M,1);


for i=1:nfd

 dL=data.map.W*(((tf+etf{i}*e-t0-et0{i}*e)*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,vdat)-...
    (tf-etf{i}*e-t0+et0{i}*e)*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,vdat)).*snm)/(2*e);
Lz=Lz+sparse(1,idx(:,i),dL,1,nz);

end



% Compute contribution of E to the gradient

Ez=spalloc(1,nz,nt+np+2*(n+m));
idx=data.FD.index.Ey;
nfd=size(idx,2);                               

et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

for i=1:nfd
Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
              E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(:,i),tf-etf(:,i),vdat))/(2*e);
end



% Return the gradient

grad=Lz+Ez;

