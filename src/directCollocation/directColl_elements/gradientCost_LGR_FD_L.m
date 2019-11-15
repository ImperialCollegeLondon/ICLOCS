function [ Lz, JL ] = gradientCost_LGR_FD_L( Lz, L, M, nz, X, Xr, U, Ur, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

