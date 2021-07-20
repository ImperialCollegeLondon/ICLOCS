function [ Lz, JL ] = gradientCost_AN_L( dL, Lz,L,nt,np,n,m,nz,X,Xr,U,Ur,P,t,T,t0,tf,vdat,data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
if np && ~isempty(dL.dp)
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
    
end

