% -------------------------------------------------------------------------
%  2. Evaluate the Jacobian of the constraints
% -------------------------------------------------------------------------

function jac=jacobianFDAdigator_LGR(f,g,avrc,X_Np1,U,P,T,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
% jacobianFDAdigator_LGR - It evaluates the Jacobian of the constraints
% with Adigator
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk



e=data.options.perturbation.J;                                 % pertubation size
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,~]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;
vdat=data.data;
X=X_Np1(1:M,:);

% Compute fgz
%------------
fgdz=sparse(const_vec_Adigator.dY_location(:,1),const_vec_Adigator.dY_location(:,2),const_vec_Adigator.dY,const_vec_Adigator.dY_size(1),const_vec_Adigator.dY_size(2));

% Compute rcz
%------------

rcz=spalloc(nrc,nz,data.map.spmatsize.jSrc);
if nrc
    [ rcz ] = jacConst_LGR_RC( rcz, nrc, nz, avrc, X_Np1, U, P, t0, tf, T, e, data );
end


% Compute bz
%------------
bz=spalloc(nb,nz,data.map.spmatsize.jSb);
if nb
    [ bz ] = jacConst_LGR_FD_B( bz, nb, nz, b, x0, xf, u0, uf, p, t0, tf, e, vdat, data );
end

% Map derivatives to the jacobian
%---------------------------------
jac=[fgdz;rcz;bz];


