function hessian=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
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
nrc=nrcl+nrcu+nrce;
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';
adjoint_g=zeros(ng,M);
adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
adjoint_g=adjoint_g';
vdat=data.data;
DT=tf-t0;
fg=vdat.functionfg;

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
%

% Compute fzz and gzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));
gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
         if j==i
              [dyn_1,gp1]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
              [dyn_2,gp2]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
              [dyn_o,go]=fg(X,U,P,DT*T+t0,vdat);
             
              fp1=(DT+etf(i)-et0(i))*dyn_1;
              fp2=(DT-etf(i)+et0(i))*dyn_2;
              fo=DT*dyn_o;
              ft=(fp2-2*fo+fp1).*adjoint_f/e2;
              gt=(gp2-2*go+gp1).*adjoint_g/e2;
         else
              [dyn_pp,gpp]=fg(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);
              [dyn_pm,gpm]=fg(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
              [dyn_mm,gmm]=fg(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
              [dyn_mp,gmp]=fg(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
             
             fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*dyn_pp;
             fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*dyn_pm;
             fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*dyn_mm;
             fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*dyn_mp;
             ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
             gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
         end
         fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft',M*n,1),nz,nz);
         g_vect=reshape(gt',M*ng,1);
         gzz=gzz+sparse(idxg(data.gAllidx,i),idxg(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
       end
    end
else
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
         if j==i
              fp1=(DT+etf(i)-et0(i))*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
              fp2=(DT-etf(i)+et0(i))*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
              fo=DT*f(X,U,P,DT*T+t0,vdat);
              ft=(fp2-2*fo+fp1).*adjoint_f/e2;
         else
            fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);   
             fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
             fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
             fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
             ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
         end
         fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
       end
    end

    if ng
        idx=data.FD.index.g;nfd=size(idx,2);
        et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
        ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;ez=data.FD.vector.g.ez;
        for i=1:nfd
           for j=1:i
            if j==i
             go=g(X,U,P,DT*T+t0,vdat);
             gp1=g(X+ex{i}*e*ez(i),U+eu{i}*e*ez(i),P+ep{i}*e*ez(i),(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
             gp2=g(X-ex{j}*e*ez(i),U-eu{j}*e*ez(i),P-ep{j}*e*ez(i),(DT-etf(j)+et0(j)).*T+t0-et0(j),vdat);
             gt=(gp2-2*go+gp1).*adjoint_g/e2;
            else
             gpp=g(X+(ex{i}+ex{j})*e*ez(i),U+(eu{i}+eu{j})*e*ez(i),P+(ep{i}+ep{j})*e*ez(i),...
                (DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);  
             gpm=g(X+(ex{i}-ex{j})*e*ez(i),U+(eu{i}-eu{j})*e*ez(i),P+(ep{i}-ep{j})*e*ez(i),...
                (DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat);
             gmm=g(X-(ex{i}+ex{j})*e*ez(i),U-(eu{i}+eu{j})*e*ez(i),P-(ep{i}+ep{j})*e*ez(i),...
                (DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
             gmp=g(X-(ex{i}-ex{j})*e*ez(i),U-(eu{i}-eu{j})*e*ez(i),P-(ep{i}-ep{j})*e*ez(i),...
                (DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat);
             gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
            end
             g_vect=reshape(gt',M*ng,1);
             gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
            end
        end
    end
end

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
% Ezz=sparse(Ezz);

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
hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
hessian=tril(hessc);


%------------- END OF CODE --------------

