function [ gz ] = jacConst_FD_G( gz, M, ng, nz, g, X, U, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.g;nfd=size(idx,2);
    et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
    ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;


    for i=1:nfd

        gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
        gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
        gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);

    end


end

