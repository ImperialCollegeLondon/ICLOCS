function [ gzz ] = hessian_CD_G( gzz, M, ng, nz, g, X, U, P, t0, T, DT, e, e2, adjoint_g, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx=data.FD.index.g;nfd=size(idx,2);
et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;%ez=data.FD.vector.g.ez;


for i=1:nfd
   for j=1:i
    if j==i
     go=g(X,U,P,DT*T+t0,vdat);
     gp1=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
     gp2=g(X-ex{j}*e,U-eu{j}*e,P-ep{j}*e,(DT-etf(j)+et0(j)).*T+t0-et0(j),vdat);
     gt=(gp2-2*go+gp1).*adjoint_g/e2;
    else
     gpp=g(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
        (DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);  
     gpm=g(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,...
        (DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat);
     gmm=g(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,...
        (DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
     gmp=g(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,...
        (DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat);
     gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
    end
     g_vect=reshape(gt',M*ng,1);
     gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
    end
end
% 
% for i=1:nfd
%    for j=1:i
%     if j==i
%      go=g(X,U,P,DT*T+t0,vdat);
%      gp1=g(X+ex{i}*e*ez(i),U+eu{i}*e*ez(i),P+ep{i}*e*ez(i),(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
%      gp2=g(X-ex{j}*e*ez(i),U-eu{j}*e*ez(i),P-ep{j}*e*ez(i),(DT-etf(j)+et0(j)).*T+t0-et0(j),vdat);
%      gt=(gp2-2*go+gp1).*adjoint_g/e2;
%     else
%      gpp=g(X+(ex{i}+ex{j})*e*ez(i),U+(eu{i}+eu{j})*e*ez(i),P+(ep{i}+ep{j})*e*ez(i),...
%         (DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);  
%      gpm=g(X+(ex{i}-ex{j})*e*ez(i),U+(eu{i}-eu{j})*e*ez(i),P+(ep{i}-ep{j})*e*ez(i),...
%         (DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat);
%      gmm=g(X-(ex{i}+ex{j})*e*ez(i),U-(eu{i}+eu{j})*e*ez(i),P-(ep{i}+ep{j})*e*ez(i),...
%         (DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
%      gmp=g(X-(ex{i}-ex{j})*e*ez(i),U-(eu{i}-eu{j})*e*ez(i),P-(ep{i}-ep{j})*e*ez(i),...
%         (DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat);
%      gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
%     end
%      g_vect=reshape(gt',M*ng,1);
%      gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
%     end
% end


end

