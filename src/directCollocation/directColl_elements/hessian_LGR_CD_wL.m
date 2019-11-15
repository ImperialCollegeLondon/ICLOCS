function [ Lzz ] = hessian_LGR_CD_wL( Lzz, nz, L, X, Xr, U, Ur, P, k0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.FD.FcnTypes.Ltype
    idx=data.FD.index.Ly;
    nfd=size(idx,2);                               

    et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
    ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;



    for i=1:nfd
    dt0_1=e*et0{i};dtf_1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
        for j=1:i
          if j==i
            Lo=data.t_segment_end.*L(X,Xr,U,Ur,P,DT/2*T+k0/2,vdat);
            Lp1=(data.t_segment_end+dtf_1/2-dt0_1/2).*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf_1-dt0_1)/2*T+(k0+dtf_1+dt0_1)/2,vdat);
            Lp2=(data.t_segment_end-dtf_1/2+dt0_1/2).*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf_1+dt0_1)/2*T+(k0-dtf_1-dt0_1)/2,vdat);
            Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
          else
            dt0_2=e*et0{j};dtf_2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
            Lpp=(data.t_segment_end+dtf_1/2-dt0_1/2+dtf_2/2-dt0_2/2).*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,vdat);
            Lpm=(data.t_segment_end+dtf_1/2-dt0_1/2-dtf_2/2+dt0_2/2).*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,vdat);
            Lmp=(data.t_segment_end-dtf_1/2+dt0_1/2+dtf_2/2-dt0_2/2).*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,vdat);
            Lmm=(data.t_segment_end-dtf_1/2+dt0_1/2-dtf_2/2+dt0_2/2).*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,vdat);
            Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
          end
           Lzz=Lzz+sparse(idx(:,i),idx(:,j),Lt,nz,nz);
        end
    end
end


end

