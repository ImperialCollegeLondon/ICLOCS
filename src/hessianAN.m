
function [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)
                               
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
t=(tf-data.t0)*T+data.k0;





[HL,HE,Hf,Hg,Hb]=hessianLagrangian(X,U,P,t,E,x0,xf,u0,uf,p,tf,data);



% Define some useful variables
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});
nz=nt+np+M*n+N*m;                           % Length of the primal variable


Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';

vdat=data.data;
t0=data.t0;
DT=tf-t0;
Tj=kron(T,ones(1,n));
% Compute fzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));  % Allocate some memory


if ~isempty(Hf)
  idx=data.FD.index.f; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  if nt
   ft=(2*df{1}.*Tj+DT*Hf{1,1}.*(Tj.^2)).*adjoint_f;  
   fzz=fzz+sparse(idx(:,1),idx(:,1),reshape(ft',M*n,1),nz,nz);  
  end 
  for i=1+nt:nfd
   for j=1:i
      if ((j==1)&&nt)
        ft=(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
      else 
        ft=DT*Hf{j,i}.*adjoint_f;
      end
      fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
    end
  end
 else
% If Hf is empty the hessian related to the dynamical system is computed numerically    
   idx=data.FD.index.f;nfd=size(idx,2);
   etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
   ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
   fo=f(X,U,P,DT*T+data.k0,vdat);
   for i=1:nfd
    fp1=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);
    for j=1:i
       if j==i;fp2=fp1;else
         fp2=f(X+ex{j}*e,U+eu{j}*e,P+ep{j}*e,(DT+etf(j)).*T+data.k0,vdat);
       end
    fpp=f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)).*T+data.k0,vdat);
    fte=DT*(fpp-fp2+fo-fp1)+etf(i)*(fpp-fp1)+etf(j)*(fpp-fp2);
    ft=fte.*adjoint_f/e2;
    fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
    end
   end
end



% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

% If there are path constraints (ng=1), their hessian is computed
% numerically whenever Hg is empty, otherwise the evaluation of the analytic espression 
% is exploited

if ng
  adjoint_g=reshape(lambda(n*M+1:n*M+ng*M)',ng,M)';
  idx=data.FD.index.g; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  if ~isempty(Hg)
  Tj=kron(T,ones(1,ng));
    if nt
      gt=(Hg{1,1}).*(Tj.^2).*adjoint_g;  
      gzz=gzz+sparse(idx(:,1),idx(:,1),reshape(gt',M*ng,1),nz,nz);
    end    
    for i=1+nt:nfd
    for j=1:i
     if ((j==1)&&nt)
      gt=(Hg{j,i}).*(Tj).*adjoint_g;  
     else   
      gt=Hg{j,i}.*adjoint_g;
     end
     gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt',M*ng,1),nz,nz);
    end
   end
  else    
   idx=data.FD.index.g;nfd=size(idx,2);
   etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
   ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
   go=g(X,U,P,DT*T+data.k0,vdat);
   for i=1:nfd
    gp1=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);

    for j=1:i
    
    if j==i;gp2=gp1;else
    gp2=g(X+ex{j}*e,U+eu{j}*e,P+ep{j}*e,(DT+etf(j)).*T+data.k0,vdat);
    end

    gpp=g(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
        (DT+etf(i)+etf(j)).*T+data.k0,vdat);
    
    gt=(gpp-gp2+go-gp1).*adjoint_g/e2; 
    gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt',M*ng,1),nz,nz);
    end
   end

  end
end







% Compute (w'L)zz
% ----------------

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
idx=data.FD.index.Ly; 
nfd=size(idx,2);  
if (~isempty(HL))
 if nt&&nfd
   Lt=(2*dL{1}.*T+DT*HL{1,1}.*(T.^2)).*data.map.w;  
   Lzz=Lzz+sparse(idx(:,1),idx(:,1),reshape(Lt',M,1),nz,nz);  
  end    
  for i=1+nt:nfd
   for j=1:i
      if ((j==1)&&nt)
        Lt=(dL{i}+DT*HL{j,i}.*T).*data.map.w;
      else   
        Lt=DT*HL{j,i}.*data.map.w;
      end
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
    end
  end
else
% If HL is empty the hessian of the cost is computed numerically
                         
etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;
Lo=DT*L(X,Xr,U,Ur,P,DT*T+data.k0,vdat);
for i=1:nfd
dt1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
Lp1=(DT+dt1)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dt1)*T+data.k0,vdat);
   for j=1:i
    dt2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};

    if j==i;Lp2=Lp1;else
        Lp2=(DT+dt2)*L(X+dx2,Xr,U+du2,Ur,P+dp2,(DT+dt2)*T+data.k0,vdat);
    end

    Lpp=(DT+dt1+dt2)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dt1+dt2)*T+data.k0,vdat);
    Lt=(Lpp-Lp2+Lo-Lp1).*data.map.w/e2; 
    Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);

    end
end

end





% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));




if ~isempty(HE)
  idx=data.FD.index.Ey; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  
for i=1:nfd
for j=1:i
  Ezz(idx(i),idx(j))=HE{j,i};
end
end

else    
idx=data.FD.index.Ey;nfd=size(idx,2);                               
etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;
Eo=E(x0,xf,u0,uf,p,tf,vdat);
for i=1:nfd
    Ep1=E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),tf+etf(:,i),vdat);

    for j=1:i
    
    if j==i;Ep2=Ep1;else
    Ep2=E(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),tf+etf(:,j),vdat);
    end

    Epp=E(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+etf(:,i)+etf(:,j),vdat);
    
    Ezz(idx(i),idx(j))=(Epp+Eo-Ep1-Ep2)/e2; 
     
    end
end

end



% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));


if nb
  adjoint=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';  
  if ~isempty(Hb)
   idx=data.FD.index.b; 
   nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
   for i=1:nfd
   for j=1:i 
      bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint',nz,nz);
   end
   end
  else   
  idx=data.FD.index.b;nfd=size(idx,2);
  etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
  ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
  exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;
  bo=b(x0,xf,u0,uf,p,tf,vdat);
  for i=1:nfd
    bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),tf+etf(:,i),vdat);
    for j=1:i
     if j==i;bp2=bp1;else
     bp2=b(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),tf+etf(:,j),vdat);
     end
     bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+etf(:,i)+etf(:,j),vdat);
     bt=(bpp-bp2+bo-bp1).*adjoint'/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz);
    end
  end

 end
end





