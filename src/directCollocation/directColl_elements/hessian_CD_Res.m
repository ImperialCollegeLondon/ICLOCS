function [ ResNormz ] = hessian_CD_Res( M, nz, nt, n, m, np, ng_eq, X, U, P, t0, T, DT, e, e2, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.FD.FcnTypes.Ltype

    dataNLP=data.resmin.dataNLP;
    data=data.resmin;
       % Compute Reszz
    % ----------------                             
    et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
    ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;%ez=dataNLP.FD.vector.Ly.ez;
    
    if nt
        ResNormz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
%         Resz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
        idx=dataNLP.FD.index.Ly;
        nfd=size(idx,2);                  
        i_st=1;
        i_end=nfd;
    else
        ResNormz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
%         Resz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
        idx=dataNLP.FD.index.Ly;
        i_st=nt+1;
        i_end=nt+np+m+n;
        idx=idx-nt;
    end
    n_res=n+ng_eq;
%     adjoint_Res=lambda(ngActive+nrc+nb+2:ngActive+nrc+nb+1+n_res);
    
    for k1=1:size(data.idx_perturb_hes,2)
        for k2=1:k1
            for i=i_st:i_end
                dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
                if k1==k2
                    j_end=i;
                else
                    j_end=i_end;
                end
                
                for j=i_st:j_end
                  if j==i && k1==k2
                    [ResNorm_o]=costResidualMin_ModeMinRes(X,U,P,DT*T+t0,data);
                    [ResNorm_p]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1),U+du1.*data.idx_perturb_hes(:,k1),P+dp1,(DT+dtf1-dt01)*T+t0+dt01,data);
                    [ResNorm_m]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1),U-du1.*data.idx_perturb_hes(:,k1),P-dp1,(DT-dtf1+dt01)*T+t0-dt01,data);
                    ResNormt=(ResNorm_m-2*ResNorm_o+ResNorm_p)/e2;
%                     Rest=(Res_m-2*Res_o+Res_p).*adjoint_Res/e2;
                  else
                    dt02=e*et0{j};dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
                    [ResNorm_pp]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1)+dx2.*data.idx_perturb_hes(:,k2),U+du1.*data.idx_perturb_hes(:,k1)+du2.*data.idx_perturb_hes(:,k2),P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,data);
                    [ResNorm_pm]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1)-dx2.*data.idx_perturb_hes(:,k2),U+du1.*data.idx_perturb_hes(:,k1)-du2.*data.idx_perturb_hes(:,k2),P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,data);
                    [ResNorm_mp]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1)+dx2.*data.idx_perturb_hes(:,k2),U-du1.*data.idx_perturb_hes(:,k1)+du2.*data.idx_perturb_hes(:,k2),P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,data);
                    [ResNorm_mm]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1)-dx2.*data.idx_perturb_hes(:,k2),U-du1.*data.idx_perturb_hes(:,k1)-du2.*data.idx_perturb_hes(:,k2),P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,data);
                    ResNormt=(ResNorm_pp-ResNorm_mp+ResNorm_mm-ResNorm_pm)/e2/4;
%                     Rest=(Res_pp-Res_mp+Res_mm-Res_pm).*adjoint_Res/e2/4;
                  end
                  
                  if (nt && i>(np+nt) && j>(np+nt)) || (~nt)
                      if (k1==1 && k2==1) || (k1==4 && k2==1)
                             if mod(data.nps,3)
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                             else
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                                if k1~=k2
                                    ResNormt(end)=[];
%                                     Rest(:,end)=[];
                                    idx2(end)=[];
                                end
                             end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif k1==4 && k2==4
                            ResNormt=ResNormt(1:3:end);
%                             Rest=Rest(:,1:3:end);
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1==2 && k2==2) || (k1==4 && k2==2) || (k1==5 && k2==2)
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            switch mod(data.nps,3)
                                 case {0,2}
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                    if (k1==5 && k2==2)
                                        ResNormt(end)=[];
%                                         Rest(:,end)=[];
                                        idx2(end)=[];
                                    end
                            end
                            ResNormz=ResNormz+sparse(max(idx1,idx2),min(idx1,idx2),ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(max(idx1,idx2),n_res,1),repelem(min(idx1,idx2),n_res,1),Rest(:),nz,nz);
    
                      elseif k1==5 && k2==5
                            ResNormt=ResNormt(2:3:end);
%                             Rest=Rest(:,2:3:end);
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1==3 && k2==3) || (k1==5 && k2==3) || (k1==6 && k2==3)
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            switch mod(data.nps,3)
                                 case {0,1}
                                    ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                    if (k1==6 && k2==3)
                                        ResNormt(end)=[];
%                                         Rest(:,end)=[];
                                        idx2(end)=[];
                                    end
                            end
                            ResNormz=ResNormz+sparse(max(idx1,idx2),min(idx1,idx2),ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(max(idx1,idx2),n_res,1),repelem(min(idx1,idx2),n_res,1),Rest(:),nz,nz);
    
    
                      elseif k1==6 && k2==6
                            ResNormt=ResNormt(3:3:end);
%                             Rest=Rest(:,3:3:end);
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1==2 && k2==1) 
                             idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end-1)+ResNormt(4:3:end)+ResNormt(5:3:end)];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end-1)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                                    idx2(end)=[];
                                 case 1
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+[ResNormt(5:3:end) 0]];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+[Rest(:,5:3:end) zeros(n_res,1)]];
                                 case 2
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+ResNormt(5:3:end)];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                             end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
    
                      elseif (k1==3 && k2==2) 
                             idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end-1)+ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,1:3:end-1)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                    idx2(end)=[];
                                 case 2
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                             end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
    
                      elseif (k1==3 && k2==1) 
                             idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
%                                     Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                                    idx2(1)=[];
                                 case 1
                                    ResNormt=ResNormt(2:3:end-2)+ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,2:3:end-2)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                    idx2(1)=[];
                                 case 2
                                    ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                    idx2(1)=[];
                                    idx1(end)=[];
                             end
                             ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
    
                      elseif (k1==6 && k2==1) 
                             idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
%                                     Rest=Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                                    idx2(1)=[];
                                 case {1,2}
                                    ResNormt=ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,3:3:end)+Rest(:,4:3:end);
                                    idx2(1)=[];
                             end
                             ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
    
                      end
                  elseif (nt && j<=(np+nt) && k2==1) %k1=1 no pertub of X and U when i<=(np+nt)
                      if i<=(np+nt) 
                          if k1==1
                              ResNormt=sum(ResNormt,2);
%                               Rest=sum(sum(Rest,2));
                              idx1=idx(1,j);
                              idx2=idx(1,i);
                              ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                               Resz=Resz+sparse(idx2,idx1,Rest,nz,nz);
                          end
                      else
                          if k1==1
                                 if mod(data.nps,3)
                                    ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
%                                     Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                    idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                 else
                                    ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
%                                     Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                    idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                 end
                                 idx1=idx(1,j)*ones(length(idx2),1);
                                 ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                  Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          elseif k1==4
                                ResNormt=ResNormt(1:3:end);
%                                 Rest=Rest(:,1:3:end);
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx1=idx(1,j)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          elseif k1==2
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx1=idx(1,j)*ones(length(idx2),1);
                                switch mod(data.nps,3)
                                     case {0,2}
                                        ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                         Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                     case 1
                                        ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                         Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                end
                                ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          elseif k1==5
                                ResNormt=ResNormt(2:3:end);
%                                 Rest=Rest(:,2:3:end);
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx1=idx(1,j)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          elseif k1==3
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx1=idx(1,j)*ones(length(idx2),1);
                                switch mod(data.nps,3)
                                     case {0,1}
                                        ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                         Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                     case 2
                                        ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                         Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                end
                                ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          elseif k1==6
                                ResNormt=ResNormt(3:3:end);
%                                 Rest=Rest(:,3:3:end);
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                                idx1=idx(1,j)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx2,n_res,1),repelem(idx1,n_res,1),Rest(:),nz,nz);
    
                          end
                      end
                  end
               end
            end
        end
    end


end
    


end

