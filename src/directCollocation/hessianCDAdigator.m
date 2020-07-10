function [Lzz,Ezz,fgzz,bzz]=hessianCDAdigator(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,const_vec_Adigator,data)
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


% Compute (w'L)zz
% ----------------

idx=data.FD.index.Ly;
nfd=size(idx,2);                               

et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;ez=data.FD.vector.Ly.ez;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

for i=1:nfd
dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i}*ez(i);dx1=e*ex{i}*ez(i); du1=e*eu{i}*ez(i);
for j=1:i
  if j==i
    Lo=DT*L(X,Xr,U,Ur,P,DT*T+t0,vdat);
    Lp1=(DT+dtf1-dt01)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf1-dt01)*T+t0+dt01,vdat);
    Lp2=(DT-dtf1+dt01)*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf1+dt01)*T+t0-dt01,vdat);
    Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
  else
    dt02=e*et0{j};dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
    Lpp=(DT+dtf1+dtf2-dt01-dt02)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,vdat);
    Lpm=(DT+dtf1-dtf2-dt01+dt02)*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,vdat);
    Lmp=(DT-dtf1+dtf2+dt01-dt02)*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,vdat);
    Lmm=(DT-dtf1-dtf2+dt01+dt02)*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,vdat);
    Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
  end
   Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
end
end



% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
idx=data.FD.index.Ey;nfd=size(idx,2);                               
et0=data.FD.vector.Ey.et0;etf=data.FD.vector.Ey.etf;ep=data.FD.vector.Ey.ep;
ex0=data.FD.vector.Ey.ex0;eu0=data.FD.vector.Ey.eu0;
exf=data.FD.vector.Ey.exf;euf=data.FD.vector.Ey.euf;ez=e*data.FD.vector.Ey.ez;

for i=1:nfd
   for j=1:i
    if j==i;
     Ep1=E(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i),tf+etf(:,i),vdat);
     Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
     Ep2=E(x0-ex0(:,i)*ez(i),xf-exf(:,i)*ez(i),u0-eu0(:,i)*ez(i),uf-euf(:,i)*ez(i),p-ep(:,i)*ez(i),t0-et0(:,i),tf-etf(:,i),vdat);
     Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
    else
     Epp=E(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
     Epm=E(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);
      
     Emp=E(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
     Emm=E(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);
      
    Ezz(idx(i),idx(j))=(Epp+Emm-Epm-Emp)/e2/4;  
    end
  end
end

% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
if nb


idx=data.FD.index.b;nfd=size(idx,2);
et0=data.FD.vector.b.et0;etf=data.FD.vector.b.etf;ep=data.FD.vector.b.ep;
ex0=data.FD.vector.b.ex0;eu0=data.FD.vector.b.eu0;
exf=data.FD.vector.b.exf;euf=data.FD.vector.b.euf;ez=e*data.FD.vector.b.ez;

adjoint=data.lambda(n*M+ngActive+nrc+(~~nb):n*M+ngActive+nrc+nb).';



for i=1:nfd
   for j=1:i
    if j==i
     bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
     bp1=b(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i)*ez(i),tf+etf(:,i)*ez(i),vdat);
     bp2=b(x0-ex0(:,j)*ez(i),xf-exf(:,j)*ez(i),u0-eu0(:,j)*ez(i),uf-euf(:,j)*ez(i),p-ep(:,j)*ez(i),t0-et0(:,i)*ez(i),tf-etf(:,j)*ez(i),vdat);
     bt=(bp2-2*bo+bp1).*adjoint'/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    else
    bpp=b(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)+exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
          uf+euf(:,i)*ez(i)+euf(:,j)*ez(j),p+ep(:,i)*ez(i)+ep(:,j)*ez(j),t0+et0(:,i)*ez(i)+et0(:,j)*ez(j),tf+etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
    bpm=b(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)-exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
          uf+euf(:,i)*ez(i)-euf(:,j)*ez(j),p+ep(:,i)*ez(i)-ep(:,j)*ez(j),t0+et0(:,i)*ez(i)-et0(:,j)*ez(j),tf+etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
    bmp=b(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)+exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
          uf-euf(:,i)*ez(i)+euf(:,j)*ez(j),p-ep(:,i)*ez(i)+ep(:,j)*ez(j),t0-et0(:,i)*ez(i)+et0(:,j)*ez(j),tf-etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
    bmm=b(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)-exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
          uf-euf(:,i)*ez(i)-euf(:,j)*ez(j),p-ep(:,i)*ez(i)-ep(:,j)*ez(j),t0-et0(:,i)*ez(i)-et0(:,j)*ez(j),tf-etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
    bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 
   
    bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    end
   end
end

end



% Return the Hessian of the Lagrangian
% -------------------------------------


% hessc=data.sigma*(Lzz+Ezz)+fgzz+bzz;
% hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
% hessian=tril(hessc);


%------------- END OF CODE --------------

