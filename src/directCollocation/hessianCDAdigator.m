function [varargout]=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
%  It evaluates the Hessian of the Lagrangian with Adigator
%
% Syntax:  hessian=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
%
% Inputs:
%    see directCollocation.m
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
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
[ng_eq,ng_neq]=deal(data.sizes{15:16});
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_fg=lambda(1:n*M+ngActive);
nrc=nrcl+nrcu+nrce;

vdat=data.data;
% t0=data.t0;
DT=tf-t0;



e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
%
% Compute fgzz
% ------------
adjoint_fg_vec=adjoint_fg(const_vec_Adigator.dYdY_location(:,1));
fgzz=sparse(const_vec_Adigator.dYdY_location(:,2),const_vec_Adigator.dYdY_location(:,3),const_vec_Adigator.dYdY.*adjoint_fg_vec,const_vec_Adigator.dYdY_size(2),const_vec_Adigator.dYdY_size(3));
% ResNormz=sparse(ResNorm_vec.dYdY_location(:,1),ResNorm_vec.dYdY_location(:,2),ResNorm_vec.dYdY,ResNorm_vec.dYdY_size(1),ResNorm_vec.dYdY_size(2));


% if nargout==2 || nargout>=5 
    % Compute (w'L)zz
    % ----------------
    if data.FD.FcnTypes.Ltype
        Lzz=spalloc(nz,nz,data.map.spmatsize.hSL);
        [ Lzz ] = hessian_CD_wL( Lzz, M, nz, L, X, Xr, U, Ur, P, t0, T, DT, e, e2, vdat, data );
    else
        Lzz=sparse(nz,nz);
    end


    % Compute Ezz
    % ------------
    if data.FD.FcnTypes.Etype
        Ezz=spalloc(nz,nz,data.map.spmatsize.hSE);
        [ Ezz ] = hessian_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data );
    else
        Ezz=sparse(nz,nz);
    end
% end

if nargout==4 || nargout>=5 

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

switch nargout
    case 2
        varargout{1}=Lzz;
        varargout{2}=Ezz;
    case 4
        varargout{1}=Lzz;
        varargout{2}=Ezz;
        varargout{3}=fgzz;
        varargout{4}=bzz; 
    case 5
%         [ ResNormz ] = hessian_CD_Res( M, nz, nt, n, m, np, ng_eq, X, U, P, t0, T, DT, e, e2, data );
        ResNorm_vec=data.ResNorm_vec;
        ResNormz=sparse(ResNorm_vec.dYdY_location(:,1),ResNorm_vec.dYdY_location(:,2),ResNorm_vec.dYdY,ResNorm_vec.dYdY_size(1),ResNorm_vec.dYdY_size(2));
        varargout{1}=Lzz;
        varargout{2}=Ezz;
        varargout{3}=fgzz;
        varargout{4}=bzz;
        varargout{5}=tril(ResNormz);

end



% Return the Hessian of the Lagrangian
% -------------------------------------


% hessc=data.sigma*(Lzz+Ezz)+fgzz+bzz;
% hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
% hessian=tril(hessc);


%------------- END OF CODE --------------

