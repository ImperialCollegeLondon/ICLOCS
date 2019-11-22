function hessian=hessianCD_LGR(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%  It evaluates the Hessian of the Lagrangian with finite diferences
%  considering central difference formula
%
% Syntax: hessian=hessianCD_LGR(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation_LGR.m
%
% Outputs:
%    hessian - hessian information
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

%------------- BEGIN CODE ---------------


% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
Gamma_1toN=reshape(data.lambda(1:M*n),M,n);
adjoint_f=Gamma_1toN;
adjoint_g=zeros(M,ng);
adjoint_g(logical(data.gActiveIdx))=data.lambda(n*M+1:n*M+ngActive);
adjoint_b=data.lambda((M*n+ngActive+nrc+1):end);

vdat=data.data;
fg=vdat.functionfg;

k0=tf+t0;
DT=tf-t0;
DTLP=repmat(data.t_segment_end,1,n);
DT_ratio_diff=repmat(data.tau_segment_ratio_diff,1,n);


e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size

% Compute fzz and gzz
% ------------



if ng && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);
    gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
    [ fzz,gzz ] =hessian_LGR_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data );
else
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);
   [ fzz ] = hessian_LGR_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data );

    if ng
        gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
        [ gzz ] = hessian_LGR_CD_G( gzz, M, ng, nz, g, X, U, P, T, k0, DT, e, e2, adjoint_g, vdat, data );
    else
        gzz=sparse(nz,nz);
    end
end
  
% Compute (w'L)zz
% ----------------
if data.FD.FcnTypes.Ltype
    Lzz=spalloc(nz,nz,data.map.spmatsize.hSL);
    [ Lzz ] = hessian_LGR_CD_wL( Lzz, nz, L, X, Xr, U, Ur, P, k0, T, DT, e, e2, vdat, data );
else
    Lzz=sparse(nz,nz);
end

% Compute Ezz
% ------------
if data.FD.FcnTypes.Etype
    Ezz=spalloc(nz,nz,data.map.spmatsize.hSE);
    [ Ezz ] = hessian_LGR_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data );
else
    Ezz=sparse(nz,nz);
end

% Compute bzz
% ------------
if nb
    bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
    [ bzz ] = hessian_LGR_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data );
else
    bzz=sparse(nz,nz);
end


% Return the Hessian of the Lagrangian
% -------------------------------------
hessian=data.sigma*(Lzz+Ezz)-fzz+gzz+bzz;


%------------- END OF CODE --------------

