function [ fz,gz,Jf ] = jacConst_LGR_FD_FG( fz, gz, M, n, ng, nz, fg, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    
    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;et=data.FD.vector.f.et;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    
    for i=1:nfd
        [dyn_p,gp]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
        [dyn_m,gm]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
        if data.options.adaptseg==1 
            fp=(repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*dyn_p;
            fm=(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*dyn_m;
            fz=fz+sparse(1:M*n,idxf(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            Jf{i}=(fp-fm)/(2*e);
        else
            fp=(DTLP+etf(i)/2*DT_ratio_diff-et0(i)/2*DT_ratio_diff).*dyn_p;
            fm=(DTLP-etf(i)/2*DT_ratio_diff+et0(i)/2*DT_ratio_diff).*dyn_m;
            ft=sparse(1:M*n,idxf(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            Jf{i}=(fp-fm)/(2*e);
            fz=fz+ft;
        end
        gz=gz+sparse(1:M*ng,idxg(:,i),reshape((gp-gm)/(2*e),M*ng,1),M*ng,nz);
    end


end

