function [ fz,Jf ] = jacConst_LGR_FD_F( fz, M, n, nz, f, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;et=data.FD.vector.f.et;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    if data.options.adaptseg==1 
        for i=1:nfd
            fp=(repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            fz=fz+sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            Jf{i}=(fp-fm)/(2*e);
        end
    else
        for i=1:nfd
            fp=(DTLP+etf(i)/2*DT_ratio_diff-et0(i)/2*DT_ratio_diff).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(DTLP-etf(i)/2*DT_ratio_diff+et0(i)/2*DT_ratio_diff).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            ft=sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            fz=fz+ft;
            Jf{i}=(fp-fm)/(2*e);
        end
    end

    
    


end

