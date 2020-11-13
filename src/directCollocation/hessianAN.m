
function [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
                               
% hessianAN - Return the Hessian of the Lagrangian when the
%             analytic option has been selected
%
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    Lzz - hessian of wL wrt z
%    Ezz - hessian of E wrt z
%    fzz - hessian of f wrt z
%    gzz - hessian of g wrt z
%    bzz - hessian of b wrt z
%
% Other m-files required: hessianLagrangian.m
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


%------------- BEGIN CODE ---------------


% Define some useful variables

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
t=(tf-t0)*T+t0;
hessianLagrangian=data.analyticDeriv.hessianLagrangian;
[HL,HE,Hf,Hg,Hb]=hessianLagrangian(X,U,P,t,E,x0,xf,u0,uf,p,t0,tf,data);

% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;                           % Length of the primal variable


Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';

vdat=data.data;
DT=tf-t0;
Tj=kron(T,ones(1,n));
fg=vdat.functionfg;

% Compute fzz and gzz
% ------------



if ng
    adjoint_g=zeros(ng,M);
    adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
    adjoint_g=adjoint_g';
end

if ~isempty(Hf) && ng && ~isempty(Hg)
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);  % Allocate some memory
    gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
    [ fzz ] = hessian_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, t0, T, e, e2, DT, adjoint_f, vdat, data );
    [ gzz ] = hessian_AN_G( gzz, Hg, M, nt, ng, nz, T, adjoint_g, data );
elseif (isempty(Hf)) && ng && (isempty(Hg)) && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);  % Allocate some memory
    gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
    [ fzz,gzz ] = hessian_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, t0, T, DT, e, e2, vdat, data );
elseif ~isempty(Hf) 
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);  % Allocate some memory
    [ fzz ] = hessian_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, t0, T, e, e2, DT, adjoint_f, vdat, data );
    if ng
        gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
        if ~isempty(Hg) 
            [ gzz ] = hessian_AN_G( gzz, Hg, M, nt, ng, nz, T, adjoint_g, data );
        else
            
            [ gzz ] = hessian_CD_G( gzz, M, ng, nz, g, X, U, P, t0, T, DT, e, e2, adjoint_g, vdat, data );
        end
    else
        gzz=sparse(nz,nz);
    end
elseif isempty(Hf)
    fzz=spalloc(nz,nz,data.map.spmatsize.hSf);  % Allocate some memory
    [ fzz ] = hessian_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, t0, T, DT, e, e2, vdat, data );
    if ng
        gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
        if ~isempty(Hg)
            [ gzz ] = hessian_AN_G( gzz, Hg, M, nt, ng, nz, T, adjoint_g, data );
        else
            [ gzz ] = hessian_CD_G( gzz, M, ng, nz, g, X, U, P, t0, T, DT, e, e2, adjoint_g, vdat, data );
        end
    else
        gzz=sparse(nz,nz);
    end
end




% Compute (w'L)zz
% ----------------



if ~isempty(HL) 
    Lzz=spalloc(nz,nz,data.map.spmatsize.hSL);
    [ Lzz ] = hessian_AN_wL( dL, Lzz, HL, M, nt, nz, T, DT, data );
else
    % If HL is empty the hessian of the cost is computed numerically
    if data.FD.FcnTypes.Ltype
        Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
        [ Lzz ] = hessian_CD_wL( Lzz, M, nz, L, X, Xr, U, Ur, P, t0, T, DT, e, e2, vdat, data );
    else
        Lzz=sparse(nz,nz);
    end
end





% Compute Ezz
% ------------



if ~isempty(HE)
    Ezz=spalloc(nz,nz,data.map.spmatsize.hSE);
    [ Ezz ] = hessian_AN_E( Ezz, HE, nt, data );
else    
    if data.FD.FcnTypes.Etype
        Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
        [ Ezz ] = hessian_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data );
    else
        Ezz=sparse(nz,nz);
    end
end


% Compute bzz
% ------------



if nb
  adjoint_b=data.lambda(n*M+ngActive+nrc+(~~nb):n*M+ngActive+nrc+nb).';
  bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
  if ~isempty(Hb)
      [ bzz ] = hessian_AN_B( bzz, Hb, nz, nt, adjoint_b, data );
  else   
      [ bzz ] =  hessian_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data );
  end
else
    bzz=sparse(nz,nz);
end





