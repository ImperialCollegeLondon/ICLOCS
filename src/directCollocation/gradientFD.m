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
% t0=data.t0;
% k0=data.k0;
% Compute contribution of L to the gradient

Lz=zeros(1,nz);
if data.FD.FcnTypes.Ltype
    [ Lz, ~ ] = gradientCost_FD_L( Lz, L, M, nz, X, Xr, U, Ur, P, t0, tf, T, e, vdat, data );
end

% Compute contribution of E to the gradient

Ez=spalloc(1,nz,nt+np+2*(n+m));
if data.FD.FcnTypes.Etype    
  [ Ez ] = gradientCost_FD_E( Ez, E, x0, xf, u0, uf, p, t0, tf, e, vdat, data );
end




% Return the gradient

grad=Lz+Ez;

