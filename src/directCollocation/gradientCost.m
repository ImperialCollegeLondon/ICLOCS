function [grad,JL]=gradientCost(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,t0,tf,data)
%GRADIENTCOST - Generate gradient of the cost  when the analytic option has been selected
%
% Syntax:  [grad,JL]=gradientCost(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    Defined in directCollocation.m
%
% Outputs:
%    solution - Data structure containing the solution
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


e=data.options.perturbation.J;                                 % Pertubation size
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});

nx=M*n;                               % Number of unknown states
nu=N*m;                               % Number of unknown controls
nz=nt+np+nx+nu;                       % Length of the optimization variable 
vdat=data.data;                       % Variables used in the functions 
% t0=data.t0;
% k0=data.k0;
t=(tf-t0)*T+t0;
gradCost=data.analyticDeriv.gradCost;

% Compute contribution of L to the gradient
Lz=zeros(1,nz);

if data.FD.index.dL.flag==1          %  data.FD.index.dE.flag==1 when the analytic expression  of the gradient for the stage cost L is supplied  
    [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data);
    idx=data.FD.index.Ly;
  if nt==1
   Lz(1)=data.map.w'*(L(X,Xr,U,Ur,P,t,vdat)+(tf-t0)*(dL.dt.*T)); 
   JaL{nt}=dL.dt(:);
  elseif nt==2
   Lz(1)=-data.map.w'*(L(X,Xr,U,Ur,P,t,vdat)+(tf-t0)*(dL.dt.*T)); 
   Lz(2)=data.map.w'*(L(X,Xr,U,Ur,P,t,vdat)+(tf-t0)*(dL.dt.*T)); 
   JaL{1}=dL.dt(:);
   JaL{2}=dL.dt(:);
  end
 if np
 for i=nt+1:nt+np
   dLp=data.map.W*(tf-t0)*dL.dp(:,i-nt);
   Lz=Lz+sparse(1,idx(:,i),dLp,1,nz);
   JaL{i}=dL.dp(:,i-nt);
 end
 end
 for  i=1:n
    dLx=data.map.W*(tf-t0)*dL.dx(:,i);
    Lz=Lz+sparse(1,idx(:,nt+np+i),dLx,1,nz);
    JaL{i+nt+np}=dL.dx(:,i);
 end  
 for i=1:m
  dLu=data.map.W*(tf-t0)*dL.du(:,i);  
  Lz=Lz+sparse(1,idx(:,nt+np+n+i),dLu,1,nz);
  JaL{i+nt+np+n}=dL.du(:,i);
 end
JL=vertcat(JaL);


else   % Numerical evalution of the gradient of the stage cost
 idx=data.FD.index.Ly;
 nfd=size(idx,2);                               
 et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
 ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;
 

 %This vector is used to adjust the size of the vector of the derivative
 %when the stage cost is identically zero 

 snm=ones(M,1);
 for i=1:nfd
  dL=((tf+etf{i}*e-t0-et0{i}*e)*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)*T+t0+et0{i}*e,vdat)-...
  (tf-etf{i}*e-t0+et0{i}*e)*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)*T+t0-et0{i}*e,vdat))/(2*e);
  Lz=Lz+sparse(1,idx(:,i),data.map.W*dL,1,nz);
  JL{i}=sparse(dL);
 end
if ~nfd
 JL=0;  
end   
end



% Compute contribution of E to the gradient

Ez=spalloc(1,nz,nt+np+2*(n+m));


if data.FD.index.dE.flag==1
  [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data);    
  Ez(data.costStruct.E==1)=[dE.dt0 dE.dtf dE.dp dE.dx0 dE.du0 dE.dxf dE.duf];

else    
    
idx=data.FD.index.Ey;
nfd=size(idx,2);                               
et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

for i=1:nfd
Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
              E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(:,i),tf-etf(:,i),vdat))/(2*e);
end


end


% Return the gradient


grad=Lz+Ez;

