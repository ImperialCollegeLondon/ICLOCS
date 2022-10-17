function [ fzz ] = hessian_CD_F( fzz, adjoint_f, M, n, nz, f, X, U, P, t0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
    
    persistent solsave; 

    if data.FD.FcnTypes.Ftype==3 && (data.ProblemTypes.FixedTime || (~data.ProblemTypes.FixedTime && ~data.FD.FcnTypes.FTRelation))
        if data.ProblemTypes.FixedTime
            if ~isfield(solsave,'ft')
                ft_save=cell(nfd,nfd);
                for i=1:nfd
                   for j=1:i
                     if j==i
                          fp1=(DT+etf(i)-et0(i))*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
                          fp2=(DT-etf(i)+et0(i))*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
                          fo=DT*f(X,U,P,DT*T+t0,vdat);
                          ft_save{i,j}=(fp2-2*fo+fp1)/e2;
                     else
                         fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);   
                         fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
                         fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
                         fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
                         ft_save{i,j}=(fpp-fpm+fmm-fmp)/e2/4;
                     end
                   end
                end  
                solsave.ft=ft_save;
            end
            
            for i=1:nfd
               for j=1:i
                  ft=solsave.ft{i,j}.*adjoint_f;
                  fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
               end
            end  
        elseif ~data.ProblemTypes.FixedTime && ~data.FD.FcnTypes.FTRelation
            if ~isfield(solsave,'f')
                fp1_save=cell(nfd,nfd);
                fp2_save=cell(nfd,nfd);
                fpp_save=cell(nfd,nfd);
                fpm_save=cell(nfd,nfd);
                fmm_save=cell(nfd,nfd);
                fmp_save=cell(nfd,nfd);
                fo_save=f(X,U,P,DT*T+t0,vdat);
                for i=1:nfd
                   for j=1:i
                     if j==i
                          fp1_save{i,j}=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
                          fp2_save{i,j}=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i)).*T+t0-et0(i),vdat);
                     else
                         fpp_save{i,j}=f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);   
                         fpm_save{i,j}=f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat); 
                         fmm_save{i,j}=f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
                         fmp_save{i,j}=f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat); 
                     end
                   end
                end  
                solsave.f.fp1_save=fp1_save;
                solsave.f.fp2_save=fp2_save;
                solsave.f.fpp_save=fpp_save;
                solsave.f.fpm_save=fpm_save;
                solsave.f.fmm_save=fmm_save;
                solsave.f.fmp_save=fmp_save;
                solsave.f.fo_save=fo_save;
            end
            
            for i=1:nfd
               for j=1:i
                 if j==i
                      fp1=(DT+etf(i)-et0(i))*solsave.f.fp1_save{i,j};
                      fp2=(DT-etf(i)+et0(i))*solsave.f.fp2_save{i,j};
                      fo=DT*solsave.f.fo_save;
                      ft=(fp2-2*fo+fp1).*adjoint_f/e2;
                 else
                    fpp=(DT+etf(i)+etf(j)-et0(i)-et0(j))*solsave.f.fpp_save{i,j};   
                     fpm=(DT+etf(i)-etf(j)-et0(i)+et0(j))*solsave.f.fpm_save{i,j}; 
                     fmm=(DT-etf(i)-etf(j)+et0(i)+et0(j))*solsave.f.fmm_save{i,j};
                     fmp=(DT-etf(i)+etf(j)+et0(i)-et0(j))*solsave.f.fmp_save{i,j}; 
                     ft=(fpp-fpm+fmm-fmp).*adjoint_f/e2/4;
                 end
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
               end
            end    
        end
    else
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
    

end

