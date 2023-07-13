function [ ResNormz ] = hessian_LGR_CD_Res( M, nz, nt, n, m, np, ng_eq, npd, nps, X_Np1, U_Np1, P, k0, T_Np1, DT, e, e2, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if data.FD.FcnTypes.Ltype

    dataNLP=data.resmin.dataNLP;
    data=data.resmin;
       % Compute Reszz
    % ----------------                             
    et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
    ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;
    
    if data.free_time
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
        i_st=1;
        i_end=m+n+np;
    end
    n_res=n+ng_eq;
%     adjoint_Res=lambda(end-n_res+1:end);
    
    idx_state=[idx(:,1:n);idx(end,1:n)+1];
    idx_pert_end=size(data.idx_perturb_hes,1);
    tau_seg_idx_state=[dataNLP.tau_seg_idx;dataNLP.tau_seg_idx(end)+1];
    for i=i_st:i_end
        dtf_1=e*etf{i};dt0_1=e*et0{i};dp1=e*ep{i};dx1=e*[ex{i};ex{i}(end,:)]; du1=e*[eu{i};eu{i}(end,:)];
        for j=i_st:i
            for k1=1:size(data.idx_perturb_hes,2)
                if i==j
                    k_end=k1;
                else
                    k_end=size(data.idx_perturb_hes,2);
                end
                
                for k2=1:k_end
                    
                 if i<=n 
                    idx1=idx_state;
                    idx_pert_end_i=idx_pert_end;
                    purb_1_x=data.idx_perturb_hes(:,k1);
                    purb_1_u=data.idx_perturb_hes(:,k1);
                 else
                    idx1=idx;
                    idx_pert_end_i=idx_pert_end-1;
                    purb_1_x=data.idx_perturb_hes(:,k1);
                    purb_1_u=[data.idx_perturb_hes(1:idx_pert_end-1,k1);data.idx_perturb_hes(idx_pert_end-1,k1)];
                 end
                 
                 if j<=n 
                    idx2=idx_state;
                    idx_pert_end_j=idx_pert_end;
                    purb_2_x=data.idx_perturb_hes(:,k2);
                    purb_2_u=data.idx_perturb_hes(:,k2);
                 else
                    idx2=idx;
                    idx_pert_end_j=idx_pert_end-1;
                    purb_2_x=data.idx_perturb_hes(:,k2);
                    purb_2_u=[data.idx_perturb_hes(1:idx_pert_end-1,k2);data.idx_perturb_hes(idx_pert_end-1,k2)];
                 end
        
                 
                  if j==i && k1==k2
                    [ResNorm_o]=costResidualMin_ModeMinRes(X_Np1,U_Np1,P,DT/2*T_Np1+k0/2,data);
                    [ResNorm_p]=costResidualMin_ModeMinRes(X_Np1+dx1.*purb_1_x,U_Np1+du1.*purb_1_u,P+dp1,(DT+dtf_1-dt0_1)/2*T_Np1+(k0+dtf_1+dt0_1)/2,data);
                    [ResNorm_m]=costResidualMin_ModeMinRes(X_Np1-dx1.*purb_1_x,U_Np1-du1.*purb_1_u,P-dp1,(DT-dtf_1+dt0_1)/2*T_Np1+(k0-dtf_1-dt0_1)/2,data);
                    ResNormt=(ResNorm_m-2*ResNorm_o+ResNorm_p)/e2;
%                     Rest=(Res_m-2*Res_o+Res_p).*adjoint_Res/e2;
                  else
                    dtf_2=e*etf{j};dt0_2=e*et0{j};dp2=e*ep{j};dx2=e*[ex{j};ex{j}(end,:)]; du2=e*[eu{j};eu{j}(end,:)];
                    [ResNorm_pp]=costResidualMin_ModeMinRes(X_Np1+dx1.*purb_1_x+dx2.*purb_2_x,U_Np1+du1.*purb_1_u+du2.*purb_2_u,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T_Np1+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,data);
                    [ResNorm_pm]=costResidualMin_ModeMinRes(X_Np1+dx1.*purb_1_x-dx2.*purb_2_x,U_Np1+du1.*purb_1_u-du2.*purb_2_u,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T_Np1+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,data);
                    [ResNorm_mp]=costResidualMin_ModeMinRes(X_Np1-dx1.*purb_1_x+dx2.*purb_2_x,U_Np1-du1.*purb_1_u+du2.*purb_2_u,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T_Np1+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,data);
                    [ResNorm_mm]=costResidualMin_ModeMinRes(X_Np1-dx1.*purb_1_x-dx2.*purb_2_x,U_Np1-du1.*purb_1_u-du2.*purb_2_u,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T_Np1+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,data);
                    ResNormt=(ResNorm_pp-ResNorm_mp+ResNorm_mm-ResNorm_pm)/e2/4;
%                     Rest=(Res_pp-Res_mp+Res_mm-Res_pm).*adjoint_Res/e2/4;
                  end
    
                  
    
                 
                  if (data.free_time && i<=(m+n) && j<=(m+n)) || (~data.free_time)    
                      if (k1==1 && k2==1)
                             [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest-1);
                                idx_rest(1)=[];
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             if mod(data.nps,3)
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             else
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             end
                             ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt(ia),nz,nz);
%                              Restl=Rest(:,ia);
%                              Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                             
                      elseif (k1==2 && k2==2) 
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            switch mod(data.nps,3)
                                 case {0,2}
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                    if length(idx1)~=length(ResNormt)
                                        if length(idx1)~=length(idx2)
                                            idx2(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                        else
                                            ResNormt(end)=[];%Rest(:,end)=[];
                                        end
                                    end
                            end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1==3 && k2==3) 
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
    
                            switch mod(data.nps,3)
                                 case {0,1}
                                    ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                    if length(idx1)~=length(ResNormt)
                                        if length(idx1)~=length(idx2)
                                            idx2(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                        else
                                            ResNormt(end)=[];%Rest(:,end)=[];
                                        end
                                    end
                            end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1>=4 && k1<=3+(max(npd)-1) && k2==1) || (k2>=4 && k2<=3+(max(npd)-1) && k1==1) %4&1
                             [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest-1);
                                idx_rest(1)=[];
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             if mod(data.nps,3)
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             else
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
%                                 Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             end
                             if length(idx1)<length(ResNormt)
                                 ResNormt(end)=[];%Rest(:,end)=[];
                             end
                             if k2==1
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                             else
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                             end
                      elseif (k1>=4 && k1<=3+(max(npd)-1) && k2==2) || (k2>=4 && k2<=3+(max(npd)-1) && k1==2) %4&2
                             if k2==2
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)))-1);
                             else
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)))-1,tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             end
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest+1);
                                idx_rest(idx_rest>nps)=[];
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             switch mod(data.nps,3)
                                 case {0,2}
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                             end
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            if k2==2
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            else
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            end               
                      elseif (k1>=3+max(npd) && k1<=3+(max(npd)-1)*2 && k2==2) || (k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 && k1==2) %(5&2)
                             [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest-1);
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             switch mod(data.nps,3)
                                 case {0,2}
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                             end
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            if k2==2
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                            else
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                            end     
    
                      elseif (k1>=3+max(npd) && k1<=3+(max(npd)-1)*2 && k2==3) || (k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 && k1==3) %5&3
                             if k2==3
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)))-1);
                             else
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)))-1,tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             end
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest+1);
                                idx_rest(idx_rest>nps)=[];
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             switch mod(data.nps,3)
                                 case {0,1}
                                    ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                             end
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            if k2==3
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            else
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            end                    
                      elseif (k1>=3+(max(npd)-1)*2+1 && k1<=3+(max(npd)-1)*3 && k2==3) || (k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 && k1==3) %(6&3)
                             [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest-1);
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             switch mod(data.nps,3)
                                 case {0,1}
                                    ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                    if (k1==6 && k2==3)
%                                         Rest(end)=[];
                                        ResNormt(:,end)=[];
                                        idx2(end)=[];
                                    end
                             end
                             idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                             idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             if k2==3
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                             else
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                             end     
    
                      elseif k1>=4 && k1<=3+(max(npd)-1) && k2>=4 && k2<=3+(max(npd)-1) %4
                            [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                            ResNormt=ResNormt(idx_rest);
%                             Rest=Rest(:,idx_rest);
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);
                             
                      elseif k1>=3+max(npd) && k1<=3+(max(npd)-1)*2 && k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 %5
                            [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                            ResNormt=ResNormt(idx_rest);
%                             Rest=Rest(:,idx_rest);
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);
    
                      elseif k1>=3+(max(npd)-1)*2+1 && k1<=3+(max(npd)-1)*3 && k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 %6
                            [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                            ResNormt=ResNormt(idx_rest);
%                             Rest=Rest(:,idx_rest);
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            ResNormz=ResNormz+sparse(idx1(ia),idx2(ib),ResNormt,nz,nz);
%                             Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);
%     
                      elseif (k1==2 && k2==1) || (k1==1 && k2==2) 
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end-1)+ResNormt(4:3:end)+ResNormt(5:3:end)];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end-1)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                                    if k2==1 && length(idx1)~=length(idx2)
                                        idx2(end)=[];
                                    elseif k1==1 && length(idx1)~=length(idx2)
                                        idx1(end)=[];
                                    end
                                 case 1
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+[ResNormt(5:3:end) 0]];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+[Rest(:,5:3:end) zeros(n_res,1)]];
                                    if k2==1 && length(idx1)~=length(idx2)
                                        idx2(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                    elseif k1==1 && length(idx1)~=length(idx2)
                                        idx1(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                    end
                                 case 2
                                    ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+ResNormt(5:3:end)];
%                                     Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                             end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                      elseif (k1==3 && k2==2) || (k1==2 && k2==3) 
                             idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                             idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end-1)+ResNormt(2:3:end)+ResNormt(3:3:end);
%                                     Rest=Rest(:,1:3:end-1)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                    if k2==2 && length(idx1)~=length(idx2)
                                        idx2(end)=[];
                                    elseif k1==2 && length(idx1)~=length(idx2)
                                        idx1(end)=[];
                                    end
                                 case 2
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                     Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                    if k1==3 && length(idx1)~=length(idx2)
                                        idx2(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                    elseif k2==3 && length(idx1)~=length(idx2)
                                        idx1(end)=[];ResNormt(end)=[];%Rest(:,end)=[];
                                    end
                             end
                             ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
    
                      elseif (k1==3 && k2==1) || (k1==1 && k2==3) 
                             idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                             idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
%                                     Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                                    if k2==1 
                                        if (i~=j && i>n) || i<=n
                                            idx2(1)=[];
                                        end
                                    elseif k1==1 && i<=n 
                                        idx1(1)=[];
                                    end
    
                                 case 1
                                    ResNormt=ResNormt(2:3:end-2)+ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,2:3:end-2)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                    if k2==1 %&& i<=n
                                        idx2(1)=[];
                                    elseif k1==1 %&& i<=n
                                        idx1(1)=[];
                                    end
                                 case 2
                                    ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                    if k2==1 %&& i<=n
                                        idx2(1)=[];
                                        if i<=n
                                            idx1(end)=[];
                                        end
                                    elseif k1==1 %&& i<=n
                                        idx1(1)=[];
                                        if j<=n
                                            idx2(end)=[];
                                        end
                                    end 
                             end
                             ResNormz=ResNormz+sparse(max(idx1,idx2),min(idx1,idx2),ResNormt,nz,nz);
%                              Resz=Resz+sparse(repelem(max(idx1,idx2),n_res,1),repelem(min(idx1,idx2),n_res,1),Rest(:),nz,nz);
    
                      elseif (k1>=3+(max(npd)-1)*2+1 && k1<=3+(max(npd)-1)*3 && k2==1) || (k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 && k1==1) %6&1
                             if k2==1
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)))-1);
                             else
                                [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)))-1,tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                             end
                             if ~isempty(idx_rest)
                                idx_rest=union(idx_rest,idx_rest+1);
                                idx_rest(idx_rest>nps)=[];
                             end
                             ResNormt(~idx_rest)=0;
%                              Rest(:,~idx_rest)=0;
                             switch mod(data.nps,3)
                                 case 0
                                    ResNormt=ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
%                                     Rest=Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                                 case {1,2}
                                    ResNormt=ResNormt(3:3:end)+ResNormt(4:3:end);
%                                     Rest=Rest(:,3:3:end)+Rest(:,4:3:end);
                             end
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            if k2==1
%                                 Restl=Rest(:,ia);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ia),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            else
%                                 Restl=Rest(:,ib);
                                ResNormz=ResNormz+sparse(max(idx1(ia),idx2(ib)),min(idx1(ia),idx2(ib)),ResNormt(ib),nz,nz);
%                                 Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                            end
    
                      end
                  elseif (data.free_time && i>(m+n) && k1==1) %k1=1 no pertub of X and U when i<=(np+nt)
                      if j>(m+n)
                          if k2==1
                              ResNormt=sum(ResNormt);
%                               Rest=sum(Rest,2);
                              idx1=idx(1,i);
                              idx2=idx(1,j);
                              ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                               Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                          end
                      else
                          if k2==1
                                 if mod(data.nps,3)
                                    ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
%                                     Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                    idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                 else
                                    ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
%                                     Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                    idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                    if length(idx2)~=length(ResNormt)
                                        ResNormt(end)=[];%Rest(:,end)=[];
                                    end
                                 end
                                 idx1=idx1(1,i)*ones(length(idx2),1);
                                  ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                   Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                          elseif k2==2
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                switch mod(data.nps,3)
                                     case {0,2}
                                        ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
%                                         Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                     case 1
                                        ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
%                                         Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                        if length(idx2)~=length(ResNormt)
                                            ResNormt(end)=[];
%                                             Rest(:,end)=[];
                                        end
                                end
                                idx1=idx1(1,i)*ones(length(idx2),1);
                                 ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                  Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                          elseif k2==3
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                switch mod(data.nps,3)
                                     case {0,1}
                                        ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
%                                         Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                     case 2
                                        ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
%                                         Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                        if length(idx2)~=length(ResNormt)
                                           ResNormt(end)=[];
%                                            Rest(:,end)=[];
                                        end
                                end
                                idx1=idx1(1,i)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                                
                          elseif k2>=4 && k2<=3+(max(npd)-1) %4
                                ResNormt=ResNormt(1:3:end);
%                                 Rest=Rest(:,1:3:end);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                idx1=idx1(1,i)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                          elseif k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 %5
                                ResNormt=ResNormt(2:3:end);
%                                 Rest=Rest(:,2:3:end);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                idx1=idx1(1,i)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
    
                          elseif k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 %6
                                ResNormt=ResNormt(3:3:end);
%                                 Rest=Rest(:,3:3:end);
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                idx1=idx1(1,i)*ones(length(idx2),1);
                                ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
%                                 Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                            
                          end
                      end
                  end
               end
            end
        end
    end
end
    


end

