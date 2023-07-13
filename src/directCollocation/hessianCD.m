function [varargout]=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%  It evaluates the Hessian of the Lagrangian with finite diferences
%  considering central difference formula
%
% Syntax:  hessian=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    hessian - hessian sparse matrix
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
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
[ng_eq,ng_neq]=deal(data.sizes{15:16});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
vdat=data.data;
DT=tf-t0;
fg=vdat.functionfg;

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
%

if nargout==2 || nargout==5 || nargout==6
    % Compute (w'L)zz
    % ----------------
    if data.FD.FcnTypes.Ltype==0 || data.FD.FcnTypes.Ltype==2
        Lzz=sparse(nz,nz);
    else
        Lzz=spalloc(nz,nz,data.map.spmatsize.hSL);
        [ Lzz ] = hessian_CD_wL( Lzz, M, nz, L, X, Xr, U, Ur, P, t0, T, DT, e, e2, vdat, data );
    end


    % Compute Ezz
    % ------------
    if data.FD.FcnTypes.Etype
        Ezz=spalloc(nz,nz,data.map.spmatsize.hSE);
        [ Ezz ] = hessian_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data );
    else
        Ezz=sparse(nz,nz);
    end
end


if nargout==3 || nargout==5 || nargout==6
    
    lambda=data.lambda(:);
    adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';
    adjoint_g=zeros(ng,M);
    adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
    adjoint_g=adjoint_g';
    % Compute fzz and gzz
    % ------------

    if ng && size(data.FD.index.f,2)==size(data.FD.index.g,2)
        if (data.FD.FcnTypes.Ftype==0 || data.FD.FcnTypes.Ftype==2) && (data.FD.FcnTypes.Gtype==0 || data.FD.FcnTypes.Gtype==2)
            fzz=sparse(nz,nz);
            gzz=sparse(nz,nz);
        else
            fzz=spalloc(nz,nz,data.map.spmatsize.hSf);
            gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
            [ fzz,gzz ] = hessian_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, t0, T, DT, e, e2, vdat, data );
        end
    else
        if data.FD.FcnTypes.Ftype==0 || data.FD.FcnTypes.Ftype==2
            fzz=sparse(nz,nz);
        else
           fzz=spalloc(nz,nz,data.map.spmatsize.hSf);
           [ fzz ] = hessian_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, t0, T, DT, e, e2, vdat, data );
        end
        if ~ng || data.FD.FcnTypes.Ftype==0 || data.FD.FcnTypes.Ftype==2
            gzz=sparse(nz,nz);
        else
            gzz=spalloc(nz,nz,data.map.spmatsize.hSg);
            [ gzz ] = hessian_CD_G( gzz, M, ng, nz, g, X, U, P, t0, T, DT, e, e2, adjoint_g, vdat, data );
        end
    end




    % Compute bzz
    % ------------

    if nb
        bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
        adjoint_b=data.lambda(n*M+ngActive+nrc+(~~nb):n*M+ngActive+nrc+nb).';
        [ bzz ] =  hessian_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data );
    else
        bzz=sparse(nz,nz);
    end
end

% Return the Hessian of the Lagrangian
% -------------------------------------
% hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
% hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
% hessian=tril(hessc);
switch nargout
    case 2
        varargout{1}=Lzz;
        varargout{2}=Ezz;
    case 3
        varargout{1}=fzz;
        varargout{2}=gzz; 
        varargout{3}=bzz; 
    case 5
        varargout{1}=Lzz;
        varargout{2}=Ezz;
        varargout{3}=fzz;
        varargout{4}=gzz; 
        varargout{5}=bzz;
    case 6
        [ ResNormz ] = hessian_CD_Res( M, nz, nt, n, m, np, ng_eq, X, U, P, t0, T, DT, e, e2, data );
        varargout{1}=Lzz;
        varargout{2}=Ezz;
        varargout{3}=fzz;
        varargout{4}=gzz; 
        varargout{5}=bzz;
        varargout{6}=tril(ResNormz);

end

%------------- END OF CODE --------------

