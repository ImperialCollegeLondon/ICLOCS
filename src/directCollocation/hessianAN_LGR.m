
function [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
% hessianAN_LGR - Return the Hessian of the Lagrangian when the analytic option has been selected
%
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
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

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
t=(tf-t0)/2*T+(tf+t0)/2;
alpha=(1-T)/2;
beta=(1+T)/2;
hessianLagrangian=data.analyticDeriv.hessianLagrangian;



[HL,HE,Hf,Hg,Hb]=hessianLagrangian(X,U,P,t,E,x0,xf,u0,uf,p,t0,tf,data);



% Define some useful variables
[nt,np,n,m,ng,nb,M]=deal(data.sizes{1:7});
nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
Gamma_1toN=reshape(data.lambda(1:M*n),M,n);
adjoint_f=Gamma_1toN;
adjoint_g=reshape(data.lambda((M*n+1):M*(n+ng)),M,ng);
adjoint_b=data.lambda((M*(n+ng)+1):end);

vdat=data.data;
k0=tf+t0;
DT=tf-t0;
DTLP=repmat(data.t_segment_end,1,n);
alpha_j=kron(alpha,ones(1,n));
beta_j=kron(beta,ones(1,n));
DT_ratio_diff=repmat(data.tau_segment_ratio_diff,1,n);

% Compute fzz and gzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));  % Allocate some memory
gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng
    adjoint_g=zeros(ng,M);
    adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
    adjoint_g=adjoint_g';
end

if ~isempty(Hf) && ng && ~isempty(Hg)
    [ fzz ] = hessian_LGR_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, T, k0, e, DTLP, adjoint_f, alpha_j, beta_j, vdat, data );
    [~,dg,~]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
    [ gzz ] = hessian_LGR_AN_G( gzz, dg, Hg, M, nt, ng, nz, T, k0, adjoint_g, data );
elseif isempty(Hf) && ng && isempty(Hg) && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    [ fzz,gzz ] = hessian_LGR_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data );
elseif ~isempty(Hf)
    [ fzz ] = hessian_LGR_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, T, k0, e, DTLP, adjoint_f, alpha_j, beta_j, vdat, data );
    if ng
        if ~isempty(Hg)
            [~,dg,~]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [ gzz ] = hessian_LGR_AN_G( gzz, dg, Hg, M, nt, ng, nz, T, k0, adjoint_g, data );
        else
            [ gzz ] = hessian_LGR_CD_G( gzz, M, ng, nz, g, X, U, P, T, k0, DT, e, e2, adjoint_g, vdat, data );
        end
    end
elseif isempty(Hf)
    [ fzz ] = hessian_LGR_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data );
    if ng
        if ~isempty(Hg)
            [~,dg,~]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [ gzz ] = hessian_LGR_AN_G( gzz, dg, Hg, M, nt, ng, nz, T, k0, adjoint_g, data );
        else
            [ gzz ] = hessian_LGR_CD_G( gzz, M, ng, nz, g, X, U, P, T, k0, DT, e, e2, adjoint_g, vdat, data );
        end
    end
end


% Compute (w'L)zz
% ----------------

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

if (~isempty(HL))
    [ Lzz ] = hessian_LGR_AN_wL( dL, Lzz, HL, M, nt, nz, alpha_j, beta_j, data );
else
    % If HL is empty the hessian of the cost is computed numerically
    [ Lzz ] = hessian_LGR_CD_wL( Lzz, nz, L, X, Xr, U, Ur, P, k0, T, DT, e, e2, vdat, data );
end



% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));

if ~isempty(HE)
    [ Ezz ] = hessian_LGR_AN_E( Ezz, HE, nz, data );
else    
    [ Ezz ] = hessian_LGR_CD_E( Ezz, E, x0, xf, u0, uf, p, t0, tf, e, e2, vdat, data );
end


% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));

if nb
  if ~isempty(Hb)
      [ bzz ] = hessian_LGR_AN_B( bzz, Hb, nz, adjoint_b, data );
  else   
      [ bzz ] = hessian_LGR_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data );
 end
end





