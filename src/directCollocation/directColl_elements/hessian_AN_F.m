function [ fzz ] = hessian_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, t0, T, e, DT, adjoint_f, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

      idx=data.FD.index.f; 
      nfd=size(idx,2);
      Tj=kron(T,ones(1,n));
      if ~isempty(df{1})
          if nt==1
           ft=(2*df{1}.*Tj+DT*Hf{1,1}.*(Tj.^2)).*adjoint_f;  
           fzz=fzz+sparse(idx(:,1),idx(:,1),reshape(ft',M*n,1),nz,nz);  
          elseif nt==2
            ft0t0=-(2*df{1}.*Tj+DT*Hf{1,1}.*(Tj.^2)).*adjoint_f;  
            ftftf=(2*df{2}.*Tj+DT*Hf{2,2}.*(Tj.^2)).*adjoint_f;  
            ft0tf=-(2*df{1}.*Tj+DT*Hf{1,2}.*(Tj.^2)).*adjoint_f;  
           fzz=fzz+sparse(idx(:,1),idx(:,1),reshape(ft0t0',M*n,1),nz,nz)+sparse(idx(:,2),idx(:,2),reshape(ftftf',M*n,1),nz,nz)+sparse(idx(:,2),idx(:,1),reshape(ft0tf',M*n,1),nz,nz);  
          end 
          for i=1+nt:nfd
           for j=1:i
              if (nt==1&&(j==1))
                    ft=(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
              elseif (nt==2&&(j<=nt))
                if j==1
                    ft=-(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
                elseif j==2
                    ft=(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
                end
              else 
                ft=DT*Hf{j,i}.*adjoint_f;
              end
              fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
            end
          end
      else
            et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
            ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

            for i=1:nt
               for j=1:nfd
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
          for i=1+nt:nfd
           for j=1+nt:i
              if (nt==1&&(j==1))
                    ft=(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
              elseif (nt==2&&(j<=nt))
                if j==1
                    ft=-(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
                elseif j==2
                    ft=(df{i}+DT*Hf{j,i}.*Tj).*adjoint_f;
                end
              else 
                ft=DT*Hf{j,i}.*adjoint_f;
              end
              fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
            end
          end
      end
    
    


end

