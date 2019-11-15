function [ gzz ] = hessian_LGR_CD_G( gzz, M, ng, nz, g, X, U, P, T, k0, DT, e, e2, adjoint_g, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx=data.FD.index.g;nfd=size(idx,2);
et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
for i=1:nfd
   for j=1:i
    if j==i
     go=g(X,U,P,DT/2*T+k0/2,vdat);
     gp1=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
     gp2=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
     gt=(gp2-2*go+gp1).*adjoint_g/e2;
    else
     gpp=g(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
        (DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);  
     gpm=g(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,...
        (DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat);
     gmm=g(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,...
        (DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
     gmp=g(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,...
        (DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat);
     gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
    end
     g_vect=reshape(gt,M*ng,1);
     gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
    end
end


end

