function hessian=resMin_hessianCD_LGR_ModeMinCost(L,g,X_Np1,U_Np1,P,T_Np1,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment,data)
%resMin_hessianCD_LGR_ModeMinCost - Hessian computation for
%integrated residual minimization (alternating method: cost
%minimization) using finite difference 
%
% Syntax:   hessian=resMin_hessianCD_LGR_ModeMinCost(L,g,X_Np1,U_Np1,P,T_Np1,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment,data)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------
lambda=data.lambda(:);sigma=data.sigma;
dataNLP=data.dataNLP;

% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive,ng_eq,ng_neq]=deal(dataNLP.sizes{1:19});
nrc=nrcl+nrcu+nrce;
if data.free_time
    nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
else
    nz=np+(M+1)*n+M*m;                              % Length of the primal variable
end

Xr=dataNLP.references.xr;Ur=dataNLP.references.ur;
adjoint_g=zeros(M,ng);
adjoint_g(logical(dataNLP.gActiveIdx))=lambda(1:ngActive);
adjoint_b=lambda((ngActive+nrc+1):(ngActive+nrc+nb));


X=X_Np1(1:M,:);
U=U_Np1(1:M,:);
T=T_Np1(1:M,:);

vdat=dataNLP.data;
k0=tf+t0;
DT=tf-t0;
DTLP=repmat(t_segment_end,1,n);


e=dataNLP.options.perturbation.H;
e2=e*e;                     % Pertubation size

% Compute gzz
% ------------
gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng
idx=dataNLP.FD.index.g;nfd=size(idx,2);
et0=e*dataNLP.FD.vector.g.et0;etf=e*dataNLP.FD.vector.g.etf;ep=dataNLP.FD.vector.g.ep;
ex=dataNLP.FD.vector.g.ex;eu=dataNLP.FD.vector.g.eu;

if data.free_time
    i_st=1;
    i_end=nfd;
else
    i_st=1;
    i_end=m+n+np;
end
for i=i_st:i_end
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
     g_vect=reshape(gt,M*ng,1);
     gzz=gzz+sparse(idx(dataNLP.gAllidx,i),idx(dataNLP.gAllidx,j),g_vect(dataNLP.gAllidx),nz,nz);
    end
end


end

  
% Compute (w'L)zz
% ----------------
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;

% 
if data.free_time
    Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
    idx=dataNLP.FD.index.Ly;
    nfd=size(idx,2);    
    i_st=1;
    i_end=nfd;
else
    Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
    idx=dataNLP.FD.index.Ly;
    i_st=1;
    i_end=m+n+np;
end

for i=i_st:i_end
dt0_1=e*et0{i};dtf_1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
    for j=1:i
      if j==i
        Lo=t_segment_end.*L(X,Xr,U,Ur,P,DT/2*T+k0/2,vdat);
        Lp1=(t_segment_end+dtf_1/2-dt0_1/2).*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf_1-dt0_1)/2*T+(k0+dtf_1+dt0_1)/2,vdat);
        Lp2=(t_segment_end-dtf_1/2+dt0_1/2).*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf_1+dt0_1)/2*T+(k0-dtf_1-dt0_1)/2,vdat);
        Lt=(Lp2-2*Lo+Lp1).*dataNLP.map.w/e2;
      else
        dt0_2=e*et0{j};dtf_2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
        Lpp=(t_segment_end+dtf_1/2-dt0_1/2+dtf_2/2-dt0_2/2).*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,vdat);
        Lpm=(t_segment_end+dtf_1/2-dt0_1/2-dtf_2/2+dt0_2/2).*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,vdat);
        Lmp=(t_segment_end-dtf_1/2+dt0_1/2+dtf_2/2-dt0_2/2).*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,vdat);
        Lmm=(t_segment_end-dtf_1/2+dt0_1/2-dtf_2/2+dt0_2/2).*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,vdat);
        Lt=(Lpp-Lmp+Lmm-Lpm).*dataNLP.map.w/e2/4;
      end
       Lzz=Lzz+sparse(idx(:,i),idx(:,j),Lt,nz,nz);
    end
end


% Compute Ezz
% ------------
et0=e*dataNLP.FD.vector.Ey.et0;etf=e*dataNLP.FD.vector.Ey.etf;ep=e*dataNLP.FD.vector.Ey.ep;
ex0=e*dataNLP.FD.vector.Ey.ex0;eu0=e*dataNLP.FD.vector.Ey.eu0;
exf=e*dataNLP.FD.vector.Ey.exf;euf=e*dataNLP.FD.vector.Ey.euf;

if data.free_time
    Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
    idx=dataNLP.FD.index.Ey;nfd=size(idx,2);                                        
    i_st=1;
    i_end=nfd;
else
    Ezz=spalloc(nz,nz,(2*m+2*n+np)*(2*m+2*n+np));
    idx=dataNLP.FD.index.Ey;
    i_st=1;
    i_end=(m+n)*2+np;
end
for i=i_st:i_end
   for j=1:i
    if j==i
     Ep1=E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat);
     Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
     Ep2=E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),t0-et0(:,i),tf-etf(:,i),vdat);
     Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
    else
     Epp=E(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
     Epm=E(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
          uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);
      
     Emp=E(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
          uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
     Emm=E(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
          uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);
    Ezz(idx(i),idx(j))=(Epp+Emm-Epm-Emp)/e2/4;  
    end
  end
end


% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
if nb 
    idx=dataNLP.FD.index.b;nfd=size(idx,2);
    et0=e*dataNLP.FD.vector.b.et0;etf=e*dataNLP.FD.vector.b.etf;ep=e*dataNLP.FD.vector.b.ep;et=e.*dataNLP.FD.vector.b.et;
    ex0=e*dataNLP.FD.vector.b.ex0;eu0=e*dataNLP.FD.vector.b.eu0;
    exf=e*dataNLP.FD.vector.b.exf;euf=e*dataNLP.FD.vector.b.euf;


    if data.free_time
        i_st=1;
        i_end=nfd;
    else
        i_st=1;
        i_end=(m+n)*2+np;
    end
    for i=i_st:i_end
       for j=1:i
        if j==i
         bo=b(x0,xf,u0,uf,p,t0,tf,vdat,dataNLP.options,t_segment);
         bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,dataNLP.options,t_segment+et(i,:)');
         bp2=b(x0-ex0(:,j),xf-exf(:,j),u0-eu0(:,j),uf-euf(:,j),p-ep(:,j),t0-et0(:,j),tf-etf(:,j),vdat,dataNLP.options,t_segment-et(j,:)');
         bt=(bp2-2*bo+bp1).*adjoint_b/e2; 
         bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        else
        bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
              uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat,dataNLP.options,t_segment+et(i,:)'+et(j,:)');
        bpm=b(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
              uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat,dataNLP.options,t_segment+et(i,:)'-et(j,:)');
        bmp=b(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
              uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat,dataNLP.options,t_segment-et(i,:)'+et(j,:)');
        bmm=b(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
              uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat,dataNLP.options,t_segment-et(i,:)'-et(j,:)');
        bt=(bpp-bpm+bmm-bmp).*adjoint_b/e2/4; 

        bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        end
       end
    end

end

% Compute Res_zz
% ----------------
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;

Resz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
idx=dataNLP.FD.index.Ly;
nfd=size(idx,2);                  
i_st=1;
i_end=nfd;
n_res=n+ng_eq;
adjoint_Res=lambda(end-n_res+1:end);

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
                [~, Res_o]=constResidualMin_ModeMinCost(X_Np1,U_Np1,P,DT/2*T_Np1+k0/2,data,1);
                [~, Res_p]=constResidualMin_ModeMinCost(X_Np1+dx1.*purb_1_x,U_Np1+du1.*purb_1_u,P+dp1,(DT+dtf_1-dt0_1)/2*T_Np1+(k0+dtf_1+dt0_1)/2,data,1);
                [~, Res_m]=constResidualMin_ModeMinCost(X_Np1-dx1.*purb_1_x,U_Np1-du1.*purb_1_u,P-dp1,(DT-dtf_1+dt0_1)/2*T_Np1+(k0-dtf_1-dt0_1)/2,data,1);
                Rest=(Res_m-2*Res_o+Res_p).*adjoint_Res/e2;
              else
                dtf_2=e*etf{j};dt0_2=e*et0{j};dp2=e*ep{j};dx2=e*[ex{j};ex{j}(end,:)]; du2=e*[eu{j};eu{j}(end,:)];
                [~, Res_pp]=constResidualMin_ModeMinCost(X_Np1+dx1.*purb_1_x+dx2.*purb_2_x,U_Np1+du1.*purb_1_u+du2.*purb_2_u,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T_Np1+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,data,1);
                [~, Res_pm]=constResidualMin_ModeMinCost(X_Np1+dx1.*purb_1_x-dx2.*purb_2_x,U_Np1+du1.*purb_1_u-du2.*purb_2_u,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T_Np1+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,data,1);
                [~, Res_mp]=constResidualMin_ModeMinCost(X_Np1-dx1.*purb_1_x+dx2.*purb_2_x,U_Np1-du1.*purb_1_u+du2.*purb_2_u,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T_Np1+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,data,1);
                [~, Res_mm]=constResidualMin_ModeMinCost(X_Np1-dx1.*purb_1_x-dx2.*purb_2_x,U_Np1-du1.*purb_1_u-du2.*purb_2_u,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T_Np1+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,data,1);
                Rest=(Res_pp-Res_mp+Res_mm-Res_pm).*adjoint_Res/e2/4;
              end

              if (data.free_time && i<=(m+n) && j<=(m+n)) || (~data.free_time)    
                  if (k1==1 && k2==1)
                         [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                         if ~isempty(idx_rest)
                            idx_rest=union(idx_rest,idx_rest-1);
                            idx_rest(1)=[];
                         end
                         Rest(:,~idx_rest)=0;
                         if mod(data.nps,3)
                            Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         else
                            Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         end
                         Restl=Rest(:,ia);
                         Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                         
                  elseif (k1==2 && k2==2) 
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        switch mod(data.nps,3)
                             case {0,2}
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                             case 1
                                Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                if length(idx1)~=size(Rest,2)
                                    if length(idx1)~=length(idx2)
                                        idx2(end)=[];Rest(:,end)=[];
                                    else
                                        Rest(:,end)=[];
                                    end
                                end
                        end
                         Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                  elseif (k1==3 && k2==3) 
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);

                        switch mod(data.nps,3)
                             case {0,1}
                                Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 2
                                Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                if length(idx1)~=size(Rest,2)
                                    if length(idx1)~=length(idx2)
                                        idx2(end)=[];Rest(:,end)=[];
                                    else
                                        Rest(:,end)=[];
                                    end
                                end
                        end
                         Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                  elseif (k1>=4 && k1<=3+(max(npd)-1) && k2==1) || (k2>=4 && k2<=3+(max(npd)-1) && k1==1) %4&1
                         [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                         if ~isempty(idx_rest)
                            idx_rest=union(idx_rest,idx_rest-1);
                            idx_rest(1)=[];
                         end
                         Rest(:,~idx_rest)=0;
                         if mod(data.nps,3)
                            Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         else
                            Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                            idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         end
                         if length(idx1)<size(Rest,2)
                             Rest(:,end)=[];
                         end
                         if k2==1
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                         else
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
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
                         Rest(:,~idx_rest)=0;
                         switch mod(data.nps,3)
                             case {0,2}
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                             case 1
                                Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                         end
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        if k2==2
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        else
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        end               
                  elseif (k1>=3+max(npd) && k1<=3+(max(npd)-1)*2 && k2==2) || (k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 && k1==2) %(5&2)
                         [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                         if ~isempty(idx_rest)
                            idx_rest=union(idx_rest,idx_rest-1);
                         end
                         Rest(:,~idx_rest)=0;
                         switch mod(data.nps,3)
                             case {0,2}
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                             case 1
                                Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                         end
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        if k2==2
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
                        else
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Restl(:),nz,nz);
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
                         Rest(:,~idx_rest)=0;
                         switch mod(data.nps,3)
                             case {0,1}
                                Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 2
                                Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                         end
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        if k2==3
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        else
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        end                    
                  elseif (k1>=3+(max(npd)-1)*2+1 && k1<=3+(max(npd)-1)*3 && k2==3) || (k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 && k1==3) %(6&3)
                         [idx_rest, ia, ib]=intersect(tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),tau_seg_idx_state(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                         if ~isempty(idx_rest)
                            idx_rest=union(idx_rest,idx_rest-1);
                         end
                         Rest(:,~idx_rest)=0;
                         switch mod(data.nps,3)
                             case {0,1}
                                Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 2
                                Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                if (k1==6 && k2==3)
                                    Rest(end)=[];
                                    idx2(end)=[];
                                end
                         end
                         idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                         idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         if k2==3
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                         else
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                         end     

                  elseif k1>=4 && k1<=3+(max(npd)-1) && k2>=4 && k2<=3+(max(npd)-1) %4
                        [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                        Rest=Rest(:,idx_rest);
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);
                         
                  elseif k1>=3+max(npd) && k1<=3+(max(npd)-1)*2 && k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 %5
                        [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                        Rest=Rest(:,idx_rest);
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);

                  elseif k1>=3+(max(npd)-1)*2+1 && k1<=3+(max(npd)-1)*3 && k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 %6
                        [idx_rest, ia, ib]=intersect(dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1))),dataNLP.tau_seg_idx(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2))));
                        Rest=Rest(:,idx_rest);
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        Resz=Resz+sparse(repelem(idx1(ia),n_res,1),repelem(idx2(ib),n_res,1),Rest(:),nz,nz);

                  elseif (k1==2 && k2==1) || (k1==1 && k2==2) 
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end-1)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                                if k2==1 && length(idx1)~=length(idx2)
                                    idx2(end)=[];
                                elseif k1==1 && length(idx1)~=length(idx2)
                                    idx1(end)=[];
                                end
                             case 1
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+[Rest(:,5:3:end) zeros(n_res,1)]];
                                if k2==1 && length(idx1)~=length(idx2)
                                    idx2(end)=[];Rest(:,end)=[];
                                elseif k1==1 && length(idx1)~=length(idx2)
                                    idx1(end)=[];Rest(:,end)=[];
                                end
                             case 2
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                         end
                         Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                  elseif (k1==3 && k2==2) || (k1==2 && k2==3) 
                         idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                         idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 1
                                Rest=Rest(:,1:3:end-1)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                if k2==2 && length(idx1)~=length(idx2)
                                    idx2(end)=[];
                                elseif k1==2 && length(idx1)~=length(idx2)
                                    idx1(end)=[];
                                end
                             case 2
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                if k1==3 && length(idx1)~=length(idx2)
                                    idx2(end)=[];Rest(:,end)=[];
                                elseif k2==3 && length(idx1)~=length(idx2)
                                    idx1(end)=[];Rest(:,end)=[];
                                end
                         end
                         Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);


                  elseif (k1==3 && k2==1) || (k1==1 && k2==3) 
                         idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                         idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                                if k2==1 
                                    if (i~=j && i>n) || i<=n
                                        idx2(1)=[];
                                    end
                                elseif k1==1 && i<=n 
                                    idx1(1)=[];
                                end

                             case 1
                                Rest=Rest(:,2:3:end-2)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                if k2==1 %&& i<=n
                                    idx2(1)=[];
                                elseif k1==1 %&& i<=n
                                    idx1(1)=[];
                                end
                             case 2
                                Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+Rest(:,4:3:end);
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
                         Resz=Resz+sparse(repelem(max(idx1,idx2),n_res,1),repelem(min(idx1,idx2),n_res,1),Rest(:),nz,nz);

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
                         Rest(:,~idx_rest)=0;
                         switch mod(data.nps,3)
                             case 0
                                Rest=Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)];
                             case {1,2}
                                Rest=Rest(:,3:3:end)+Rest(:,4:3:end);
                         end
                        idx1=idx1(logical(data.idx_perturb_hes(1:idx_pert_end_i,k1)),i);
                        idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                        if k2==1
                            Restl=Rest(:,ia);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        else
                            Restl=Rest(:,ib);
                            Resz=Resz+sparse(repelem(max(idx1(ia),idx2(ib)),n_res,1),repelem(min(idx1(ia),idx2(ib)),n_res,1),Restl(:),nz,nz);
                        end

                  end
              elseif (data.free_time && i>(m+n) && k1==1) %k1=1 no pertub of X and U when i<=(np+nt)
                  if j>(m+n)
                      if k2==1
                          Rest=sum(Rest,2);
                          idx1=idx(1,i);
                          idx2=idx(1,j);
                          Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                      end
                  else
                      if k2==1
                             if mod(data.nps,3)
                                Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                             else
                                Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n_res,1)]];
                                idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                                if length(idx2)~=size(Rest,2)
                                    Rest(:,end)=[];
                                end
                             end
                             idx1=idx1(1,i)*ones(length(idx2),1);
                             Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                      elseif k2==2
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            switch mod(data.nps,3)
                                 case {0,2}
                                    Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n_res,1)];
                                    if length(idx2)~=size(Rest,2)
                                        Rest(:,end)=[];
                                    end
                            end
                            idx1=idx1(1,i)*ones(length(idx2),1);
                            Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                      elseif k2==3
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            switch mod(data.nps,3)
                                 case {0,1}
                                    Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n_res,1)];
                                    if length(idx2)~=size(Rest,2)
                                       Rest(:,end)=[];
                                    end
                            end
                            idx1=idx1(1,i)*ones(length(idx2),1);
                            Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);
                            
                      elseif k2>=4 && k2<=3+(max(npd)-1) %4
                            Rest=Rest(:,1:3:end);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            idx1=idx1(1,i)*ones(length(idx2),1);
                            Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                      elseif k2>=3+max(npd) && k2<=3+(max(npd)-1)*2 %5
                            Rest=Rest(:,2:3:end);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            idx1=idx1(1,i)*ones(length(idx2),1);
                            Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                      elseif k2>=3+(max(npd)-1)*2+1 && k2<=3+(max(npd)-1)*3 %6
                            Rest=Rest(:,3:3:end);
                            idx2=idx2(logical(data.idx_perturb_hes(1:idx_pert_end_j,k2)),j);
                            idx1=idx1(1,i)*ones(length(idx2),1);
                            Resz=Resz+sparse(repelem(idx1,n_res,1),repelem(idx2,n_res,1),Rest(:),nz,nz);

                      end
                  end
              end
           end
        end
    end
end

% Return the Hessian of the Lagrangian
% -------------------------------------?
hessc=gzz+bzz+Resz;
hessian=sigma*(Lzz+Ezz)+hessc;
