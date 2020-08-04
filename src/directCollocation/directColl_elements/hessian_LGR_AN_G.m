function [ gzz ] = hessian_LGR_AN_G( gzz, dg, Hg, M, nt, ng, nz, X, U, P, T, DT, k0, adjoint_g, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  idx=data.FD.index.g; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
    if dg.flag && ~isempty(dg.dt)
          for i=1:nfd
            for j=1:i
              if j==(nfd-nt+1) && length(df)==nfd
                 gt=(Hg{j,i}).*(alpha_j).*adjoint_g;   
              elseif j==(nfd) && length(df)==nfd
                 gt=(Hg{j,i}).*(beta_j).*adjoint_g;  
              else 
                 gt=Hg{j,i}.*adjoint_g;
              end
             gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt,M*ng,1),nz,nz);
            end
          end
    else
        for i=1:nfd-nt
          for j=1:i
           gt=Hg{j,i}.*adjoint_g;
           gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt,M*ng,1),nz,nz);
          end
        end
      
        idx=data.FD.index.g;nfd=size(idx,2);
        et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
        ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
        for i=nfd-nt+1:nfd
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
             gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt,M*ng,1),nz,nz);
           end
        end
    end


end

