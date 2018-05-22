
% -------------------------------------------------------------------------
%   Evaluate the Hessian of the Lagrangian
% -------------------------------------------------------------------------

function hessian=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,fxx_adigator,data)

%  It evaluates the Hessian of the Lagrangian with Adigator
%
% Syntax:  hessian=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,fxx_adigator,data)
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
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';
adjoint_g=reshape(lambda(n*M+1:n*M+ng*M)',ng,M)';
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

if nt
    for i=1:nfd
       for j=1:i
         if i~=1 && j~=1 && i<=(nfd-np-nb) && j<=(nfd-np-nb)
             if ~isempty(fxx_adigator)
                  ft=DT*fxx_adigator{i-1,j-1}'.*adjoint_f;
             else
                  ft=zeros(M*n,1);
             end
         else
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
         end
         fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
       end
    end
else
    for i=1:nfd
       for j=1:i
           if i<=(nfd-np-nb) && j<=(nfd-np-nb)
             if ~isempty(fxx_adigator)
                  ft=DT*fxx_adigator{i,j}'.*adjoint_f;
             else
                  ft=zeros(M*n,1);
             end
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
end


% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng
idx=data.FD.index.g;nfd=size(idx,2);
etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
for i=1:nfd
   for j=1:i
    if j==i;
     go=g(X,U,P,DT*T+data.k0,vdat);
     gp1=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);
     gp2=g(X-ex{j}*e,U-eu{j}*e,P-ep{j}*e,(DT-etf(j)).*T+data.k0,vdat);
     gt=(gp2-2*go+gp1).*adjoint_g/e2;
    else
     gpp=g(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
        (DT+etf(i)+etf(j)).*T+data.k0,vdat);  
     gpm=g(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,...
        (DT+etf(i)-etf(j)).*T+data.k0,vdat);
     gmm=g(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,...
        (DT-etf(i)-etf(j)).*T+data.k0,vdat);
     gmp=g(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,...
        (DT-etf(i)+etf(j)).*T+data.k0,vdat);
     gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
    end
     gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt',M*ng,1),nz,nz);
    end
end


end


% Compute (w'L)zz
% ----------------

idx=data.FD.index.Ly;
nfd=size(idx,2);                               

etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

for i=1:nfd
dt1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
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
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

for i=1:nfd
   for j=1:i
    if j==i;
     Ep1=E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0,tf+etf(:,i),vdat);
     Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
     Ep2=E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0,tf-etf(:,i),vdat);
     Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
    else
     Epp=E(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0,tf+etf(:,i)+etf(:,j),vdat);
     Epm=E(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
          uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0,tf+etf(:,i)-etf(:,j),vdat);
      
     Emp=E(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
          uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0,tf-etf(:,i)+etf(:,j),vdat);
     Emm=E(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
          uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0,tf-etf(:,i)-etf(:,j),vdat);
      
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
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;

adjoint=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';



for i=1:nfd
   for j=1:i
    if j==i;
     bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
     bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0,tf+etf(:,i),vdat);
     bp2=b(x0-ex0(:,j),xf-exf(:,j),u0-eu0(:,j),uf-euf(:,j),p-ep(:,j),t0,tf-etf(:,j),vdat);
     bt=(bp2-2*bo+bp1).*adjoint'/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    else
    bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0,tf+etf(:,i)+etf(:,j),vdat);
    bpm=b(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
          uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0,tf+etf(:,i)-etf(:,j),vdat);
    bmp=b(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
          uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0,tf-etf(:,i)+etf(:,j),vdat);
    bmm=b(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
          uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0,tf-etf(:,i)-etf(:,j),vdat);
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

