function [ Lzz ] = hessian_CD_wL( Lzz, M, nz, L, X, Xr, U, Ur, P, t0, T, DT, e, e2, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    
    idx=data.FD.index.Ly;
    nfd=size(idx,2);                               

    et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
    ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;%ez=data.FD.vector.Ly.ez;
    
    persistent solsave; 
    if isfield(vdat,'wLSave') && ~isfield(solsave,'Lt') && ~isfield(solsave,'L')
        solsave=vdat.wLSave;
    end
    if isfield(solsave,'Lt') && (size(solsave.Lt,1)~=nfd || size(solsave.Lt,2)~=nfd || size(solsave.Lt{1},1)~=M)
        solsave = rmfield(solsave,'Lt');
    end
    if isfield(solsave,'L') && (size(solsave.L.Lp1_save,1)~=nfd || size(solsave.L.Lp1_save,2)~=nfd)
        solsave = rmfield(solsave,'L');
    end

    if data.FD.FcnTypes.Ltype==3 && (data.ProblemTypes.FixedTime || (~data.ProblemTypes.FixedTime && ~data.FD.FcnTypes.FTRelation))
        if data.ProblemTypes.FixedTime
            if ~isfield(solsave,'Lt')
                Lt_save=cell(nfd,nfd);
                for i=1:nfd
                   dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
                   for j=1:i
                     if j==i
                          Lo=DT*L(X,Xr,U,Ur,P,DT*T+t0,vdat);
                          Lp1=(DT+dtf1-dt01)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf1-dt01)*T+t0+dt01,vdat);
                          Lp2=(DT-dtf1+dt01)*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf1+dt01)*T+t0-dt01,vdat);
                          Lt_save{i,j}=(Lp2-2*Lo+Lp1).*data.map.w/e2;
                     else
                         dt02=e*et0{j};dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
                         Lpp=(DT+dtf1+dtf2-dt01-dt02)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,vdat);
                         Lpm=(DT+dtf1-dtf2-dt01+dt02)*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,vdat);
                         Lmp=(DT-dtf1+dtf2+dt01-dt02)*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,vdat);
                         Lmm=(DT-dtf1-dtf2+dt01+dt02)*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,vdat);
                         Lt_save{i,j}=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
                     end
                   end
                end  
                solsave.Lt=Lt_save;
            end
            
            if isfield(data.options,'parfor') && data.options.parfor
                parfor i=1:nfd
                   for j=1:i
                      Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(solsave.Lt{i,j}',M,1),nz,nz);
                   end
                end  
            else
                for i=1:nfd
                   for j=1:i
                      Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(solsave.Lt{i,j}',M,1),nz,nz);
                   end
                end  
            end
        elseif ~data.ProblemTypes.FixedTime && ~data.FD.FcnTypes.FTRelation
            if ~isfield(solsave,'L')
                Lp1_save=cell(nfd,nfd);
                Lp2_save=cell(nfd,nfd);
                Lpp_save=cell(nfd,nfd);
                Lpm_save=cell(nfd,nfd);
                Lmm_save=cell(nfd,nfd);
                Lmp_save=cell(nfd,nfd);
                Lo_save=L(X,Xr,U,Ur,P,DT*T+t0,vdat);
                for i=1:nfd
                   dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
                   for j=1:i
                     if j==i
                         Lp1_save{i,j}=L(X+dx1,Xr,U+du1,Ur,P+dp1,T,vdat);
                         Lp2_save{i,j}=L(X-dx1,Xr,U-du1,Ur,P-dp1,T,vdat);
                     else
                         dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
                         Lpp_save{i,j}=L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,T,vdat);   
                         Lpm_save{i,j}=L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,T,vdat);
                         Lmm_save{i,j}=L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,T,vdat);
                         Lmp_save{i,j}=L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,T,vdat);
                     end
                   end
                end  
                solsave.L.Lp1_save=Lp1_save;
                solsave.L.Lp2_save=Lp2_save;
                solsave.L.Lpp_save=Lpp_save;
                solsave.L.Lpm_save=Lpm_save;
                solsave.L.Lmm_save=Lmm_save;
                solsave.L.Lmp_save=Lmp_save;
                solsave.L.Lo_save=Lo_save;
            end
            
            for i=1:nfd
               dt01=e*et0{i};dtf1=e*etf{i};
               for j=1:i
                 if j==i
                      Lp1=(DT+dtf1-dt01)*solsave.L.Lp1_save{i,j};
                      Lp2=(DT-dtf1+dt01)*solsave.L.Lp2_save{i,j};
                      Lo=DT*solsave.L.Lo_save;
                      Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
                 else
                     dt02=e*et0{j};dtf2=e*etf{j};
                     Lpp=(DT+dtf1+dtf2-dt01-dt02)*solsave.L.Lpp_save{i,j};   
                     Lpm=(DT+dtf1-dtf2-dt01+dt02)*solsave.L.Lpm_save{i,j}; 
                     Lmm=(DT-dtf1+dtf2+dt01-dt02)*solsave.L.Lmm_save{i,j};
                     Lmp=(DT-dtf1-dtf2+dt01+dt02)*solsave.L.Lmp_save{i,j}; 
                     Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
                 end
                 Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
               end
            end    
        end
    else
        for i=1:nfd
            dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
            for j=1:i
              if j==i
                Lo=DT*L(X,Xr,U,Ur,P,DT*T+t0,vdat);
                Lp1=(DT+dtf1-dt01)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf1-dt01)*T+t0+dt01,vdat);
                Lp2=(DT-dtf1+dt01)*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf1+dt01)*T+t0-dt01,vdat);
                Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
              else
                dt02=e*et0{j};dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
                Lpp=(DT+dtf1+dtf2-dt01-dt02)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,vdat);
                Lpm=(DT+dtf1-dtf2-dt01+dt02)*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,vdat);
                Lmp=(DT-dtf1+dtf2+dt01-dt02)*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,vdat);
                Lmm=(DT-dtf1-dtf2+dt01+dt02)*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,vdat);
                Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
              end
               Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
            end
        end
    end
    
    
    global wLSave;
    if isempty(wLSave)
        wLSave=solsave;
    end
    
    
    
end
    



