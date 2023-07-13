function [varargout]=hessianCDAdigator_LGR(L,f,g,X_Np1,U_Np1,P,T,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)

%  It evaluates the Hessian of the Lagrangian with Adigator
%
% Syntax:  hessian=hessianCDAdigator_LGR(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
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
% Gamma_1toN=reshape(data.lambda(1:M*n),M,n);
% D_N=diag(1./data.map.w);
% adjoint_f=D_N*Gamma_1toN;

adjoint_fg=data.lambda(1:n*M+ngActive);
% 
% adjoint_f=Gamma_1toN;
% adjoint_g=zeros(M,ng);
% adjoint_g(logical(data.gActiveIdx))=data.lambda(n*M+1:n*M+ngActive);

adjoint_b=data.lambda((M*n+ngActive+nrc+1):end);




vdat=data.data;
% t0=data.t0;
k0=tf+t0;
DT=tf-t0;
DTLP=repmat(data.t_segment_end,1,n);
% DTLP=(tf-t0)/2;
X=X_Np1(1:end-1,:);
U=U_Np1(1:end-1,:);

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
%


% Compute fgzz
% ------------
adjoint_fg_vec=adjoint_fg(const_vec_Adigator.dYdY_location(:,1));
fgzz=sparse(const_vec_Adigator.dYdY_location(:,2),const_vec_Adigator.dYdY_location(:,3),const_vec_Adigator.dYdY.*adjoint_fg_vec,const_vec_Adigator.dYdY_size(2),const_vec_Adigator.dYdY_size(3));


  
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

if nargout==4 || nargout>=5
% Compute bzz
% ------------
    if nb
        bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
        [ bzz ] = hessian_LGR_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data );
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
%         [ng_eq,ng_neq]=deal(data.resmin.dataNLP.sizes{18:19});
%         [ ResNormz ] = hessian_LGR_CD_Res( M, nz, nt, n, m, np, ng_eq, npd, nps, X_Np1, U_Np1, P, k0, [T;1], DT, e, e2, data );
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
% hessian=data.sigma*(Lzz+Ezz)+fgzz+bzz;

% spy(hessc)
% figure
% hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';



%------------- END OF CODE --------------