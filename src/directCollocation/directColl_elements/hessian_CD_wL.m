function [ Lzz ] = hessian_CD_wL( Lzz, M, nz, L, X, Xr, U, Ur, P, t0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.FD.FcnTypes.Ltype
    
    idx=data.FD.index.Ly;
    nfd=size(idx,2);                               

    et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
    ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;ez=data.FD.vector.Ly.ez;


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
end
    


end

