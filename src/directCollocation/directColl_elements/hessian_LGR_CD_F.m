function [ fzz ] = hessian_LGR_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
             if j==i
                 fp1=(DTLP-et0(i)/2*DT_ratio_diff+etf(i)/2*DT_ratio_diff).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                 fp2=(DTLP+et0(i)/2*DT_ratio_diff-etf(i)/2*DT_ratio_diff).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                 fo=DTLP.*f(X,U,P,DT/2*T+k0/2,vdat);
                 ft=(fp2-2*fo+fp1)/e2;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             else
                 fpp=(DTLP+(-et0(i)-et0(j)+etf(i)+etf(j))/2*DT_ratio_diff).*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);   
                 fpm=(DTLP+(-et0(i)+et0(j)+etf(i)-etf(j))/2*DT_ratio_diff).*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                 fmm=(DTLP+(+et0(i)+et0(j)-etf(i)-etf(j))/2*DT_ratio_diff).*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                 fmp=(DTLP+(+et0(i)-et0(j)-etf(i)+etf(j))/2*DT_ratio_diff).*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                 ft=(fpp-fpm+fmm-fmp)/e2/4;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             end
       end
    end
    


end

