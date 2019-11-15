function [ Ez ] = gradientCost_LGR_FD_E( Ez, E, x0, xf, u0, uf, p, t0, tf, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx=data.FD.index.Ey;
nfd=size(idx,2);                               
et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

for i=1:nfd
Ez(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat)-...
              E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0+et0(:,i),tf-etf(:,i),vdat))/(2*e);
end


end

