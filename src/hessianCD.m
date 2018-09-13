
% -------------------------------------------------------------------------
%   Evaluate the Hessian of the Lagrangian
% -------------------------------------------------------------------------

function hessian=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)

%  It evaluates the Hessian of the Lagrangian with finite diferences
%  considering central difference formula
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------


% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{:});
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
% adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)'
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';
adjoint_g=zeros(ng,M);
adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
adjoint_g=adjoint_g';
% adjoint_g=reshape(lambda(n*M+1:n*M+ng*M)',ng,M)';
% adjoint_rc=reshape(lambda((n+ng)*M+1:(n+ng)*M+nrc)',nrc,1)';
% adjoint_rc=lambda(n*M+M*ng+1:n*M+M*ng+nrc).';
% adjoint_b=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';
vdat=data.data;
t0=data.t0;
DT=tf-t0;



e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
%
% Compute fzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));


idx=data.FD.index.f;nfd=size(idx,2);
etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

for i=1:nfd
   for j=1:i
     if j==i
          fp1=(DT+etf(i))*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);
          fp2=(DT-etf(i))*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)).*T+data.k0,vdat);
          fo=DT*f(X,U,P,DT*T+data.k0,vdat);
          ft=(fp2-2*fo+fp1).*adjoint_f/e2;
     else
        fpp=(DT+etf(i)+etf(j))*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)).*T+data.k0,vdat);   
         fpm=(DT+etf(i)-etf(j))*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)).*T+data.k0,vdat); 
         fmm=(DT-etf(i)-etf(j))*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)).*T+data.k0,vdat);
         fmp=(DT-etf(i)+etf(j))*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)).*T+data.k0,vdat); 
         ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
     end
     fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
   end
end


% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng
idx=data.FD.index.g;nfd=size(idx,2);
etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;ez=data.FD.vector.g.ez;
for i=1:nfd
   for j=1:i
    if j==i
     go=g(X,U,P,DT*T+data.k0,vdat);
     gp1=g(X+ex{i}*e*ez(i),U+eu{i}*e*ez(i),P+ep{i}*e*ez(i),(DT+etf(i)).*T+data.k0,vdat);
     gp2=g(X-ex{j}*e*ez(i),U-eu{j}*e*ez(i),P-ep{j}*e*ez(i),(DT-etf(j)).*T+data.k0,vdat);
     gt=(gp2-2*go+gp1).*adjoint_g/e2;
    else
     gpp=g(X+(ex{i}+ex{j})*e*ez(i),U+(eu{i}+eu{j})*e*ez(i),P+(ep{i}+ep{j})*e*ez(i),...
        (DT+etf(i)+etf(j)).*T+data.k0,vdat);  
     gpm=g(X+(ex{i}-ex{j})*e*ez(i),U+(eu{i}-eu{j})*e*ez(i),P+(ep{i}-ep{j})*e*ez(i),...
        (DT+etf(i)-etf(j)).*T+data.k0,vdat);
     gmm=g(X-(ex{i}+ex{j})*e*ez(i),U-(eu{i}+eu{j})*e*ez(i),P-(ep{i}+ep{j})*e*ez(i),...
        (DT-etf(i)-etf(j)).*T+data.k0,vdat);
     gmp=g(X-(ex{i}-ex{j})*e*ez(i),U-(eu{i}-eu{j})*e*ez(i),P-(ep{i}-ep{j})*e*ez(i),...
        (DT-etf(i)+etf(j)).*T+data.k0,vdat);
     gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
    end
     g_vect=reshape(gt',M*ng,1);
     gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
    end
end


end

% Compute rczz
% ------------

% rczz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

% if nrc
% idx=data.FD.index.rc(1:M:end,:);nfd=size(idx,2);
% etf=e*data.FD.vector.rc.etf;ep=data.FD.vector.rc.ep;
% ex=data.FD.vector.rc.ex;eu=data.FD.vector.rc.eu;
% for i=1:nfd
%    for j=1:i
%     if j==i;
%      rco=avrc(X,U,P,DT*T+data.k0,vdat);
%      rcp1=avrc(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);
%      rcp2=avrc(X-ex{j}*e,U-eu{j}*e,P-ep{j}*e,(DT-etf(j)).*T+data.k0,vdat);
%      rct=(rcp2-2*rco+rcp1).*adjoint_rc/e2;
%     else
%      rcpp=avrc(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
%         (DT+etf(i)+etf(j)).*T+data.k0,vdat);  
%      rcpm=avrc(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,...
%         (DT+etf(i)-etf(j)).*T+data.k0,vdat);
%      rcmm=avrc(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,...
%         (DT-etf(i)-etf(j)).*T+data.k0,vdat);
%      rcmp=avrc(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,...
%         (DT-etf(i)+etf(j)).*T+data.k0,vdat);
%      rct=(rcpp-rcpm+rcmm-rcmp).*adjoint_rc/e2/4;
%     end
%      rczz=rczz+sparse(idx(:,i),idx(:,j),rct,nz,nz);
%     end
% end
% end

% Compute (w'L)zz
% ----------------

idx=data.FD.index.Ly;
nfd=size(idx,2);                               

etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;ez=data.FD.vector.Ly.ez;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

for i=1:nfd
dt1=e*etf{i};dp1=e*ep{i}*ez(i);dx1=e*ex{i}*ez(i); du1=e*eu{i}*ez(i);
for j=1:i
  if j==i;
    Lo=DT*L(X,Xr,U,Ur,P,DT*T+data.k0,vdat);
    Lp1=(DT+dt1)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dt1)*T+data.k0,vdat);
    Lp2=(DT-dt1)*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dt1)*T+data.k0,vdat);
    Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
  else
    dt2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
    Lpp=(DT+dt1+dt2)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dt1+dt2)*T+data.k0,vdat);
    Lpm=(DT+dt1-dt2)*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dt1-dt2)*T+data.k0,vdat);
    Lmp=(DT-dt1+dt2)*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dt1+dt2)*T+data.k0,vdat);
    Lmm=(DT-dt1-dt2)*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dt1-dt2)*T+data.k0,vdat);
    Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
  end
   Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
end
end



% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
idx=data.FD.index.Ey;nfd=size(idx,2);                               
etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;ez=e*data.FD.vector.Ey.ez;

for i=1:nfd
   for j=1:i
    if j==i;
     Ep1=E(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0,tf+etf(:,i),vdat);
     Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
     Ep2=E(x0-ex0(:,i)*ez(i),xf-exf(:,i)*ez(i),u0-eu0(:,i)*ez(i),uf-euf(:,i)*ez(i),p-ep(:,i)*ez(i),t0,tf-etf(:,i),vdat);
     Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
    else
     Epp=E(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0,tf+etf(:,i)+etf(:,j),vdat);
     Epm=E(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0,tf+etf(:,i)-etf(:,j),vdat);
      
     Emp=E(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0,tf-etf(:,i)+etf(:,j),vdat);
     Emm=E(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0,tf-etf(:,i)-etf(:,j),vdat);
      
    Ezz(idx(i),idx(j))=(Epp+Emm-Epm-Emp)/e2/4;  
    end
  end
end






% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
if nb


idx=data.FD.index.b;nfd=size(idx,2);
etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;ez=e*data.FD.vector.b.ez;

adjoint=data.lambda(n*M+ngActive+nrc+(~~nb):n*M+ngActive+nrc+nb).';



for i=1:nfd
   for j=1:i
    if j==i;
     bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
     bp1=b(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0,tf+etf(:,i)*ez(i),vdat);
     bp2=b(x0-ex0(:,j)*ez(i),xf-exf(:,j)*ez(i),u0-eu0(:,j)*ez(i),uf-euf(:,j)*ez(i),p-ep(:,j)*ez(i),t0,tf-etf(:,j)*ez(i),vdat);
     bt=(bp2-2*bo+bp1).*adjoint'/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    else
    bpp=b(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0,tf+etf(:,i)+etf(:,j),vdat);
    bpm=b(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0,tf+etf(:,i)-etf(:,j),vdat);
    bmp=b(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0,tf-etf(:,i)+etf(:,j),vdat);
    bmm=b(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0,tf-etf(:,i)-etf(:,j),vdat);
    bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 
   
    bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    end
   end
end

end


% Return the Hessian of the Lagrangian
% -------------------------------------


hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
hessian=tril(hessc);


%------------- END OF CODE --------------

