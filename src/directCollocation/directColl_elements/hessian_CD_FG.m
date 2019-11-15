function [ fzz,gzz ] = hessian_CD_FG( fzz, gzz, adjoint_f, adjoint_g, M, n, ng, nz, fg, X, U, P, t0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
         if j==i
              [dyn_1,gp1]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
              [dyn_2,gp2]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
              [dyn_o,go]=fg(X,U,P,DT*T+t0,vdat);
             
              fp1=(DT+etf(i)-et0(i))*dyn_1;
              fp2=(DT-etf(i)+et0(i))*dyn_2;
              fo=DT*dyn_o;
              ft=(fp2-2*fo+fp1).*adjoint_f/e2;
              gt=(gp2-2*go+gp1).*adjoint_g/e2;
         else
              [dyn_pp,gpp]=fg(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);
              [dyn_pm,gpm]=fg(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
              [dyn_mm,gmm]=fg(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
              [dyn_mp,gmp]=fg(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
             
             fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*dyn_pp;
             fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*dyn_pm;
             fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*dyn_mm;
             fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*dyn_mp;
             ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
             gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
         end
         fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft',M*n,1),nz,nz);
         g_vect=reshape(gt',M*ng,1);
         gzz=gzz+sparse(idxg(data.gAllidx,i),idxg(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
       end
    end


end

