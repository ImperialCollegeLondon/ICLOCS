function [ fzz ] = hessian_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, t0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
         if j==i
              fp1=(DT+etf(i)-et0(i))*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
              fp2=(DT-etf(i)+et0(i))*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
              fo=DT*f(X,U,P,DT*T+t0,vdat);
              ft=(fp2-2*fo+fp1).*adjoint_f/e2;
         else
            fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);   
             fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
             fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
             fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
             ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
         end
         fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
       end
    end    
    


end

