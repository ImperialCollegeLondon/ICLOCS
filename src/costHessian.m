function  HC=costHessian(z,data,dL)

% COSTHESSIAN - Return the Hessian of the cost for fmincon
%
% Syntax:  [HC] =costHessian(z,data) 
%
% Inputs:
%    z     - Unknown NLP vector
%    data  - Data passed to the functions evaluated during optimization
%
% Outputs:
%    HC    - Hessian of the cost at z
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

%------------- BEGIN CODE --------------





% Define some useful variables

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size



% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{:});
nz=nt+np+M*n+N*m;                           % Length of the primal variable
[L,E,f,g,b]=deal(data.functions{:});
vdat=data.data;

% Get matrices
mp=data.map;

%--------------------------------------------------------------------------
% Extract and format vectors from NLP variable
%--------------------------------------------------------------------------

% Extract states and inputs from z and reshape for function evaluations
X=reshape(mp.Vx*z,n,M)';
usp=reshape(mp.Vu*z,m,N)';
U=kron(usp,ones((M-1)/N,1));
U=[U;U(end,:)];



% Extract design parameters if specified and convert to cells
if np; P=reshape(repmat(z(nt+1:nt+np),M,1),np,M)';else P=spalloc(M,0,1);end
% check this line

% Construct time vector
if nt; tf=z(1); else tf=data.tf; end
t0=data.t0;
k0=data.k0;

% if strcmp(data.options.transcription,'discrete'); tf=1; t0=0; end
T=[0;cumsum(data.tau)]*data.Nm/ns;

% Extract x0,u0,xf,uf,p
p=z(nt+1:nt+np);
x0=z(nt+np+1:nt+np+n);
u0=z(nt+np+(M-1)/N*n+1:nt+np+(M-1)/N*n+m);
xf=z(end-n+1:end);
uf=z(end-m-n+1:end-n);

% Format reference inputs and states if applicable
Xr=data.references.xr;Ur=data.references.ur;


vdat=data.data;
DT=tf-t0;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));

if strcmp(data.options.derivatives,'analytic')
  [HL,HE]=hessianLagrangian(X,U,P,(tf-t0)*T+k0,E,x0,xf,u0,uf,p,tf,data);
  % Compute Lzz
  % ------------
  if ~isempty(HL)
   idx=data.FD.index.Ly; 
   nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
   if nt
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
     if j<i
      ltz=reshape(Lt',M,1);
      Lzz=Lzz+sparse(idx(:,i),idx(:,j),ltz,nz,nz);
      Lzz=Lzz+sparse(idx(:,j),idx(:,i),ltz,nz,nz);
     else    
      Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
     end
    end
   end
 else
 % If HL is empty the hessian of the cost is computed numerically
  idx=data.FD.index.Ly;
  nfd=size(idx,2);                               
  etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
  ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;
  Lo=DT*L(X,Xr,U,Ur,P,DT*T+data.k0,vdat);
  for i=1:nfd
    dt1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
    Lp1=(DT+dt1)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dt1)*T+k0,vdat);
   for j=1:i
    dt2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};

    if j==i;Lp2=Lp1;else
        Lp2=(DT+dt2)*L(X+dx2,Xr,U+du2,Ur,P+dp2,(DT+dt2)*T+k0,vdat);
    end

    Lpp=(DT+dt1+dt2)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dt1+dt2)*T+k0,vdat);
    Lt=(Lpp-Lp2+Lo-Lp1).*data.map.w/e2;
    if j<i
     ltz=reshape(Lt',M,1);
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),ltz,nz,nz);
     Lzz=Lzz+sparse(idx(:,j),idx(:,i),ltz,nz,nz);
    else
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);   
    end
   end
  end
  end
  % Compute Ezz
  % ------------

  if ~isempty(HE)
   idx=data.FD.index.Ey; 
   nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
   for i=1:nfd
   for j=1:i
      Ezz(idx(i),idx(j))=HE{j,i};
    if j<i
      Ezz(idx(j),idx(i))=HE{j,i};  
    end
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
    
    if j<i
      Ezzij=(Epp+Eo-Ep1-Ep2)/e2;   
      Ezz(idx(j),idx(i))=Ezzij;  
      Ezz(idx(i),idx(j))=Ezzij; 
    else
      Ezz(idx(i),idx(j))=(Epp+Eo-Ep1-Ep2)/e2;  
    end
    end
   end
 end
else  %if strcmp(data.options.derivatives,'numeric')
  %%%%%%
  % Compute Lzz (numeric option)
  %%%%%
  idx=data.FD.index.Ly;
  nfd=size(idx,2);                               
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
    if j<i
     ltz=reshape(Lt',M,1);
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),ltz,nz,nz);
     Lzz=Lzz+sparse(idx(:,j),idx(:,i),ltz,nz,nz);
    else
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);   
    end
    
    
   end
  end
  
  %%%%%%
  % Compute Ezz (numeric option)
  %%%%%
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
    
    if j<i
      Ezzij=(Epp+Eo-Ep1-Ep2)/e2;   
      Ezz(idx(j),idx(i))=Ezzij;  
      Ezz(idx(i),idx(j))=Ezzij; 
    else
      Ezz(idx(i),idx(j))=(Epp+Eo-Ep1-Ep2)/e2;  
    end
 
    end
    end 
  
  
  
end


HC=Lzz+Ezz;