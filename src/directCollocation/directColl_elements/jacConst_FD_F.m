function [ fz,Jf ] = jacConst_FD_F( fz, M, n, nz, f, X, U, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;


    for i=1:nfd
       fp=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
       fm=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
       ft=((tf+etf(i)-t0-et0(i))*fp-(tf-etf(i)-t0+et0(i))*fm)/(2*e);
       fz=fz+sparse(1:M*n,idx(:,i),reshape(ft',M*n,1),M*n,nz);
       Jf{i}=(fp-fm)/(2*e);
    end
    
    


end

