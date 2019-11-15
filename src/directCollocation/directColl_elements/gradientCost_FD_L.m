function [ Lz, JL ] = gradientCost_FD_L( Lz, L, M, nz, X, Xr, U, Ur, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

