function [ rcz ] = jacConst_RC( rcz, nrc, nz, avrc, X, U, P, t0, tf, T, e, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



    idx=data.FD.index.rc;nfd=size(idx,2);
    et0=e*data.FD.vector.rc.et0;etf=e*data.FD.vector.rc.etf;ep=data.FD.vector.rc.ep;
    ex=data.FD.vector.rc.ex;eu=data.FD.vector.rc.eu;

    for i=1:nfd
        if ~any(ex{i}(:)) && ~any(eu{i}(:))
            rcp=avrc(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),data);
            rcm=avrc(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),data);
            rcz=rcz+sparse(1:nrc,idx(:,i),(rcp-rcm)/(2*e),nrc,nz);
        end
    end
    rcz=rcz+[data.map.Acl;data.map.Ae;data.map.Acu];




end

