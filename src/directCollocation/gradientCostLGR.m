function [grad,JL]=gradientCostLGR(L,X,Xr,U,Ur,P,T,t,E,x0,xf,u0,uf,p,t0,tf,data)
%gradientCostLGR - Generate gradient of the cost  when the analytic option has been selected
%
% Syntax:  [grad,JL]=gradientCostLGR(L,X,Xr,U,Ur,P,T,t,E,x0,xf,u0,uf,p,t0,tf,data)
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
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});

nx=(M+1)*n;                               % Number of unknown states
nu=M*m;                               % Number of unknown controls
nz=nt+np+nx+nu;                       % Length of the optimization variable 
vdat=data.data;                       % Variables used in the functions 
gradCost=data.analyticDeriv.gradCost;

% t0=data.t0;
% k0=data.k0;

% Compute contribution of L to the gradient

Lz=zeros(1,nz);

if data.FD.index.dL.flag==1          %  data.FD.index.dE.flag==1 when the analytic expression  of the gradient for the stage cost L is supplied  
    [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data);
    
    
    
    idx=data.FD.index.Ly;
    for  i=1:n
       dLx=data.map.w.*(data.t_segment_end.*dL.dx(:,i));
       Lz=Lz+sparse(1,idx(:,i),dLx,1,nz);
       JaL{i}=dL.dx(:,i);
    end  
    

    for i=1:m
       dLu=data.map.w.*(data.t_segment_end.*dL.du(:,i));  
       Lz=Lz+sparse(1,idx(:,n+i),dLu,1,nz);
       JaL{i+n}=dL.du(:,i);
    end
    

    if np
        for i=1:np
            dLp=data.map.w*(data.t_segment_end.*dL.dp(:,i));
            Lz=Lz+sparse(1,idx(:,n+m+i),dLp,1,nz);
            JaL{i+n+m}=dL.dp(:,i);
        end
    end
    

    if nt==2 && ~isempty(dL.dt)
        for i=1:nt
            if i==1 
                dLt=-0.5*data.map.w'*L(X,Xr,U,Ur,P,t,vdat)+data.map.w'*(data.t_segment_end.*(dL.dt.*(1-T)/2)); %t0
            elseif i==2 
                dLt=0.5*data.map.w'*L(X,Xr,U,Ur,P,t,vdat)+data.map.w'*(data.t_segment_end.*(dL.dt.*(1+T)/2)); %tf
            end
            Lz=Lz+sparse(1,idx(:,n+m+np+i),dLt,1,nz);
            JaL{i+n+m+np}=dL.dt(:);
        end
    end
   
    
    JL=vertcat(JaL);

else   % Numerical evalution of the gradient of the stage cost
 idx=data.FD.index.Ly;
 nfd=size(idx,2);   
 et=data.FD.vector.Ly.et;et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;
 ep=data.FD.vector.Ly.ep;
 ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;
 

 %This vector is used to adjust the size of the vector of the derivative
 %when the stage cost is identically zero 

 if data.options.adaptseg==1 
     snm=ones(M,1);
     for i=1:nfd
      dL=(((data.t_segment_end+data.t_segment_mat_m*(et{i}'*e)).*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*T+(tf+etf{i}*e+t0+et0{i}*e)/2,vdat)-...
      (data.t_segment_end-data.t_segment_mat_m*(et{i}'*e)).*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*T+(tf-etf{i}*e+t0-et0{i}*e)/2,vdat)).*snm)/(2*e);
      Lz=Lz+sparse(1,idx(:,i),diag(data.map.w)*dL,1,nz);
      JL{i}=dL;
     end
 else
     snm=ones(M,1);
     for i=1:nfd
      dL=(((data.t_segment_end+etf{i}*e/2-et0{i}*e/2).*L(X+ex{i}*e,Xr,U+eu{i}*e,Ur,P+ep{i}*e,(tf+etf{i}*e-t0-et0{i}*e)/2*T+(tf+etf{i}*e+t0+et0{i}*e)/2,vdat)-...
      (data.t_segment_end-etf{i}*e/2+et0{i}*e/2).*L(X-ex{i}*e,Xr,U-eu{i}*e,Ur,P-ep{i}*e,(tf-etf{i}*e-t0+et0{i}*e)/2*T+(tf-etf{i}*e+t0-et0{i}*e)/2,vdat)).*snm)/(2*e);
      Lz=Lz+sparse(1,idx(:,i),data.map.w.*dL,1,nz);
      JL{i}=dL;
     end
 end
 
if ~nfd
 JL=0;  
end   
end



% Compute contribution of E to the gradient

Ez=spalloc(1,nz,nt+np+2*(n+m));


if data.FD.index.dE.flag==1
    [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data);    
    Ez(data.costStruct.E==1)=[dE.dx0 dE.dxf dE.du0 dE.duf dE.dp dE.dt0 dE.dtf ];

else    
    
    idx=data.FD.index.Ey;
    nfd=size(idx,2);                               
    et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
    ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
    exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

    for i=1:nfd
    Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
                  E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat))/(2*e);
    end

   
end





% Return the gradient



grad=Lz+Ez;


