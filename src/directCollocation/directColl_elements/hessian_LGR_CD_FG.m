function [ fzz,gzz ] = hessian_LGR_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, T, k0, DTLP, DT, DT_ratio_diff, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
    
    for i=1:nfd
       for j=1:i
             if j==i
                 [dyn_1,gp1]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                 [dyn_2,gp2]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                 [dyn_o,go]=fg(X,U,P,DT/2*T+k0/2,vdat);
                 
                 fp1=(DTLP-et0(i)/2*DT_ratio_diff+etf(i)/2*DT_ratio_diff).*dyn_1;
                 fp2=(DTLP+et0(i)/2*DT_ratio_diff-etf(i)/2*DT_ratio_diff).*dyn_2;
                 fo=DTLP.*dyn_o;
                 ft=(fp2-2*fo+fp1)/e2;
                 fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 gt=(gp2-2*go+gp1).*adjoint_g/e2;
             else
                 [dyn_pp,gpp]=fg(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);
                 [dyn_pm,gpm]=fg(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                 [dyn_mm,gmm]=fg(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                 [dyn_mp,gmp]=fg(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                 
                 fpp=(DTLP+(-et0(i)-et0(j)+etf(i)+etf(j))/2*DT_ratio_diff).*dyn_pp;   
                 fpm=(DTLP+(-et0(i)+et0(j)+etf(i)-etf(j))/2*DT_ratio_diff).*dyn_pm; 
                 fmm=(DTLP+(+et0(i)+et0(j)-etf(i)-etf(j))/2*DT_ratio_diff).*dyn_mm;
                 fmp=(DTLP+(+et0(i)-et0(j)-etf(i)+etf(j))/2*DT_ratio_diff).*dyn_mp; 
                 ft=(fpp-fpm+fmm-fmp)/e2/4;
                 fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
             end
             g_vect=reshape(gt,M*ng,1);
             gzz=gzz+sparse(idxg(data.gAllidx,i),idxg(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
       end
    end

end

