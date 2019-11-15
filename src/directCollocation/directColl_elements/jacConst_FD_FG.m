function [ fz,gz,Jf ] = jacConst_FD_FG( fz, gz, M, n, ng, nz, fg, X, U, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
    
    for i=1:nfd
        [dyn_p,gp]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
        [dyn_m,gm]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
        fp=(tf+etf(i)-t0-et0(i))*dyn_p;
        fm=(tf-etf(i)-t0+et0(i))*dyn_m;
        Jfi=(fp-fm)'/(2*e);
        fz=fz+sparse(1:M*n,idxf(:,i),reshape(Jfi,M*n,1),M*n,nz);
        gz=gz+sparse(1:M*ng,idxg(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);
        Jf{i}=Jfi';
    end


end

