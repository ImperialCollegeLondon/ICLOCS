function hessian=resMin_hessianCD_ModeMinWeightedCostRes(L,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%resMin_hessianCD_LGR_ModeMinWeightedCostRes - Hessian computation for
%integrated residual minimization (weighted cost method) using finite difference
%
% Syntax:   hessian=resMin_hessianCD_LGR_ModeMinWeightedCostRes(L,g,X_Np1,U_Np1,P,T_Np1,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment,data)
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
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(dataNLP.sizes{1:13});
nrc=nrcl+nrcu+nrce;
if nt
    nz=M*n+N*m+np+nt;                           % Length of the primal variable
else
    nz=M*n+N*m+np;                           % Length of the primal variable
end
Xr=dataNLP.references.xr;Ur=dataNLP.references.ur;
adjoint_g=zeros(ng,M);
adjoint_g(logical(dataNLP.gActiveIdx'))=lambda(1:ngActive);
adjoint_g=adjoint_g';
vdat=dataNLP.data;
DT=tf-t0;



e=dataNLP.options.perturbation.H;
e2=e*e;                     % Pertubation size


% Compute gzz
% ------------
if nt
    gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));
else
    gzz=spalloc(nz,nz,M*((m+n)*(m+n)+np));
end

if ng
    if nt
        idx=dataNLP.FD.index.g;
        nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        idx=dataNLP.FD.index.g;
        i_st=nt+1;
        i_end=nt+np+m+n;
        idx=idx-nt;
    end
    et0=e*dataNLP.FD.vector.g.et0;etf=e*dataNLP.FD.vector.g.etf;ep=dataNLP.FD.vector.g.ep;
    ex=dataNLP.FD.vector.g.ex;eu=dataNLP.FD.vector.g.eu;ez=dataNLP.FD.vector.g.ez;

    for i=i_st:i_end
       for j=i_st:i
        if j==i
         go=g(X,U,P,DT*T+t0,vdat);
         gp1=g(X+ex{i}*e*ez(i),U+eu{i}*e*ez(i),P+ep{i}*e*ez(i),(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
         gp2=g(X-ex{j}*e*ez(i),U-eu{j}*e*ez(i),P-ep{j}*e*ez(i),(DT-etf(j)+et0(j)).*T+t0-et0(i),vdat);
         gt=(gp2-2*go+gp1).*adjoint_g/e2;
        else
         gpp=g(X+(ex{i}+ex{j})*e*ez(i),U+(eu{i}+eu{j})*e*ez(i),P+(ep{i}+ep{j})*e*ez(i),...
            (DT+etf(i)+etf(j)-et0(i)-et0(j)).*T+t0+et0(i)+et0(j),vdat);  
         gpm=g(X+(ex{i}-ex{j})*e*ez(i),U+(eu{i}-eu{j})*e*ez(i),P+(ep{i}-ep{j})*e*ez(i),...
            (DT+etf(i)-etf(j)-et0(i)+et0(j)).*T+t0+et0(i)-et0(j),vdat);
         gmm=g(X-(ex{i}+ex{j})*e*ez(i),U-(eu{i}+eu{j})*e*ez(i),P-(ep{i}+ep{j})*e*ez(i),...
            (DT-etf(i)-etf(j)+et0(i)+et0(j)).*T+t0-et0(i)-et0(j),vdat);
         gmp=g(X-(ex{i}-ex{j})*e*ez(i),U-(eu{i}-eu{j})*e*ez(i),P-(ep{i}-ep{j})*e*ez(i),...
            (DT-etf(i)+etf(j)+et0(i)-et0(j)).*T+t0-et0(i)+et0(j),vdat);
         gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
        end
         g_vect=reshape(gt',M*ng,1);
         gzz=gzz+sparse(idx(dataNLP.gAllidx,i),idx(dataNLP.gAllidx,j),g_vect(dataNLP.gAllidx),nz,nz);
        end
    end
end

% Compute (w'L)zz
% ----------------
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;ez=dataNLP.FD.vector.Ly.ez;

if nt
    Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
    idx=dataNLP.FD.index.Ly;
    nfd=size(idx,2);                  
    i_st=1;
    i_end=nfd;
else
    Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
    idx=dataNLP.FD.index.Ly;
    i_st=nt+1;
    i_end=nt+np+m+n;
    idx=idx-nt;
end
for i=i_st:i_end
    dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i}*ez(i);dx1=e*ex{i}*ez(i); du1=e*eu{i}*ez(i);
    for j=i_st:i
      if j==i;
        Lo=DT*L(X,Xr,U,Ur,P,DT*T+t0,vdat);
        Lp1=(DT+dtf1)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf1-dt01)*T+t0+dt01,vdat);
        Lp2=(DT-dtf1)*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf1+dt01)*T+t0-dt01,vdat);
        Lt=(Lp2-2*Lo+Lp1).*dataNLP.map.w/e2;
      else
        dt02=e*et0{j};dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
        Lpp=(DT+dtf1+dtf2)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,vdat);
        Lpm=(DT+dtf1-dtf2)*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,vdat);
        Lmp=(DT-dtf1+dtf2)*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,vdat);
        Lmm=(DT-dtf1-dtf2)*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,vdat);
        Lt=(Lpp-Lmp+Lmm-Lpm).*dataNLP.map.w/e2/4;
      end
       Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
    end
end



% Compute Ezz
% ------------
et0=dataNLP.FD.vector.Ey.et0;etf=dataNLP.FD.vector.Ey.etf;ep=dataNLP.FD.vector.Ey.ep;
ex0=dataNLP.FD.vector.Ey.ex0;eu0=dataNLP.FD.vector.Ey.eu0;
exf=dataNLP.FD.vector.Ey.exf;euf=dataNLP.FD.vector.Ey.euf;ez=e*dataNLP.FD.vector.Ey.ez;

if nt
    Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
    idx=dataNLP.FD.index.Ey;nfd=size(idx,2);                               
    i_st=1;
    i_end=nfd;
else
    Ezz=spalloc(nz,nz,(2*m+2*n+np)*(2*m+2*n+np));
    idx=dataNLP.FD.index.Ey;
    i_st=nt+1;
    i_end=nt+np+(m+n)*2;
    idx=idx-nt;
end
for i=i_st:i_end
   for j=i_st:i
    if j==i;
     Ep1=E(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i),tf+etf(:,i),vdat);
     Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
     Ep2=E(x0-ex0(:,i)*ez(i),xf-exf(:,i)*ez(i),u0-eu0(:,i)*ez(i),uf-euf(:,i)*ez(i),p-ep(:,i)*ez(i),t0-et0(:,i),tf-etf(:,i),vdat);
     Ezz(idx(i),idx(j))=(Ep1-2*Eo+Ep2)/e2; 
    else
     Epp=E(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
     Epm=E(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);
      
     Emp=E(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
     Emm=E(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
          uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);
      
    Ezz(idx(i),idx(j))=(Epp+Emm-Epm-Emp)/e2/4;  
    end
  end
end


% Compute bzz
% ------------
if nt
    bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
else
    bzz=spalloc(nz,nz,(2*m+2*n+np)*(2*m+2*n+np));
end

if nb && (~isfield(dataNLP.options,'resminRep') || ~dataNLP.options.resminRep.collmatch)
    if nt
        idx=dataNLP.FD.index.b;nfd=size(idx,2);
        i_st=1;
        i_end=nfd;
    else
        idx=dataNLP.FD.index.b;
        i_st=nt+1;
        i_end=nt+np+(m+n)*2;
        idx=idx-nt;
    end
    et0=dataNLP.FD.vector.b.et0;etf=dataNLP.FD.vector.b.etf;ep=dataNLP.FD.vector.b.ep;
    ex0=dataNLP.FD.vector.b.ex0;eu0=dataNLP.FD.vector.b.eu0;
    exf=dataNLP.FD.vector.b.exf;euf=dataNLP.FD.vector.b.euf;ez=e*dataNLP.FD.vector.b.ez;

    adjoint=lambda(ngActive+nrc+(~~nb):ngActive+nrc+nb).';

    for i=i_st:i_end
       for j=i_st:i
        if j==i
         bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
         bp1=b(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i)*ez(i),tf+etf(:,i)*ez(i),vdat);
         bp2=b(x0-ex0(:,j)*ez(i),xf-exf(:,j)*ez(i),u0-eu0(:,j)*ez(i),uf-euf(:,j)*ez(i),p-ep(:,j)*ez(i),t0-et0(:,i)*ez(i),tf-etf(:,j)*ez(i),vdat);
         bt=(bp2-2*bo+bp1).*adjoint'/e2; 
         bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        else
        bpp=b(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)+exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
              uf+euf(:,i)*ez(i)+euf(:,j)*ez(i),p+ep(:,i)*ez(i)+ep(:,j)*ez(i),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
        bpm=b(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf+exf(:,i)*ez(i)-exf(:,j)*ez(i),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
              uf+euf(:,i)*ez(i)-euf(:,j)*ez(i),p+ep(:,i)*ez(i)-ep(:,j)*ez(i),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);
        bmp=b(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)+exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(i),...
              uf-euf(:,i)*ez(i)+euf(:,j)*ez(i),p-ep(:,i)*ez(i)+ep(:,j)*ez(i),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
        bmm=b(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(i),xf-exf(:,i)*ez(i)-exf(:,j)*ez(i),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(i),...
              uf-euf(:,i)*ez(i)-euf(:,j)*ez(i),p-ep(:,i)*ez(i)-ep(:,j)*ez(i),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);
        bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 

        bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        end
       end
    end

end

% Compute Reszz
% ----------------                             
et0=dataNLP.FD.vector.Ly.et0;etf=dataNLP.FD.vector.Ly.etf;ep=dataNLP.FD.vector.Ly.ep;
ex=dataNLP.FD.vector.Ly.ex;eu=dataNLP.FD.vector.Ly.eu;ez=dataNLP.FD.vector.Ly.ez;

if nt
    ResNormz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
    Resz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
    idx=dataNLP.FD.index.Ly;
    nfd=size(idx,2);                  
    i_st=1;
    i_end=nfd;
else
    ResNormz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
    Resz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+np*np);
    idx=dataNLP.FD.index.Ly;
    i_st=nt+1;
    i_end=nt+np+m+n;
    idx=idx-nt;
end
adjoint_Res=lambda(ngActive+nrc+nb+2:ngActive+nrc+nb+1+n);

for k1=1:size(data.idx_perturb_hes,2)
    for k2=1:k1
        for i=i_st:i_end
            dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i}*ez(i);dx1=e*ex{i}*ez(i); du1=e*eu{i}*ez(i);
            if k1==k2
                j_end=i;
            else
                j_end=i_end;
            end
            
            for j=i_st:j_end
              if j==i && k1==k2
                [ResNorm_o, Res_o]=costResidualMin_ModeMinRes(X,U,P,DT*T+t0,data);
                [ResNorm_p, Res_p]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1),U+du1.*data.idx_perturb_hes(:,k1),P+dp1,(DT+dtf1-dt01)*T+t0+dt01,data);
                [ResNorm_m, Res_m]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1),U-du1.*data.idx_perturb_hes(:,k1),P-dp1,(DT-dtf1+dt01)*T+t0-dt01,data);
                ResNormt=(ResNorm_m-2*ResNorm_o+ResNorm_p)/e2;
                Rest=(Res_m-2*Res_o+Res_p).*adjoint_Res/e2;
              else
                dtf2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
                [ResNorm_pp, Res_pp]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1)+dx2.*data.idx_perturb_hes(:,k2),U+du1.*data.idx_perturb_hes(:,k1)+du2.*data.idx_perturb_hes(:,k2),P+dp1+dp2,(DT+dtf1+dtf2-dt01-dt02)*T+t0+dt01+dt02,data);
                [ResNorm_pm, Res_pm]=costResidualMin_ModeMinRes(X+dx1.*data.idx_perturb_hes(:,k1)-dx2.*data.idx_perturb_hes(:,k2),U+du1.*data.idx_perturb_hes(:,k1)-du2.*data.idx_perturb_hes(:,k2),P+dp1-dp2,(DT+dtf1-dtf2-dt01+dt02)*T+t0+dt01-dt02,data);
                [ResNorm_mp, Res_mp]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1)+dx2.*data.idx_perturb_hes(:,k2),U-du1.*data.idx_perturb_hes(:,k1)+du2.*data.idx_perturb_hes(:,k2),P-dp1+dp2,(DT-dtf1+dtf2+dt01-dt02)*T+t0-dt01+dt02,data);
                [ResNorm_mm, Res_mm]=costResidualMin_ModeMinRes(X-dx1.*data.idx_perturb_hes(:,k1)-dx2.*data.idx_perturb_hes(:,k2),U-du1.*data.idx_perturb_hes(:,k1)-du2.*data.idx_perturb_hes(:,k2),P-dp1-dp2,(DT-dtf1-dtf2+dt01+dt02)*T+t0-dt01-dt02,data);
                ResNormt=(ResNorm_pp-ResNorm_mp+ResNorm_mm-ResNorm_pm)/e2/4;
                Rest=(Res_pp-Res_mp+Res_mm-Res_pm).*adjoint_Res/e2/4;
              end
              
              if (nt && i>(np+nt) && j>(np+nt)) || (~nt)
                  if (k1==1 && k2==1) || (k1==4 && k2==1)
                         if mod(data.nps,3)
                            ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
                            Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                         else
                            ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
                            Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n,1)]];
                            idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                            if k1~=k2
                                ResNormt(end)=[];
                                Rest(:,end)=[];
                                idx2(end)=[];
                            end
                         end
                         ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                         Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);

                  elseif k1==4 && k2==4
                        ResNormt=ResNormt(1:3:end);
                        Rest=Rest(:,1:3:end);
                        idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                        idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                        ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                        Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);

                  elseif (k1==2 && k2==2) || (k1==4 && k2==2) || (k1==5 && k2==2)
                        idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                        idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                        switch mod(data.nps,3)
                             case {0,2}
                                ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                             case 1
                                ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
                                Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n,1)];
                                if (k1==5 && k2==2)
                                    ResNormt(end)=[];
                                    Rest(:,end)=[];
                                    idx2(end)=[];
                                end
                        end
                        ResNormz=ResNormz+sparse(max(idx1,idx2),min(idx1,idx2),ResNormt,nz,nz);
                        Resz=Resz+sparse(repelem(max(idx1,idx2),n,1),repelem(min(idx1,idx2),n,1),Rest(:),nz,nz);

                  elseif k1==5 && k2==5
                        ResNormt=ResNormt(2:3:end);
                        Rest=Rest(:,2:3:end);
                        idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                        idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                        ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                        Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);

                  elseif (k1==3 && k2==3) || (k1==5 && k2==3) || (k1==6 && k2==3)
                        idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                        idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                        switch mod(data.nps,3)
                             case {0,1}
                                ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
                                Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 2
                                ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
                                Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n,1)];
                                if (k1==6 && k2==3)
                                    ResNormt(end)=[];
                                    Rest(:,end)=[];
                                    idx2(end)=[];
                                end
                        end
                        ResNormz=ResNormz+sparse(max(idx1,idx2),min(idx1,idx2),ResNormt,nz,nz);
                        Resz=Resz+sparse(repelem(max(idx1,idx2),n,1),repelem(min(idx1,idx2),n,1),Rest(:),nz,nz);


                  elseif k1==6 && k2==6
                        ResNormt=ResNormt(3:3:end);
                        Rest=Rest(:,3:3:end);
                        idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                        idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                        ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                        Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);

                  elseif (k1==2 && k2==1) 
                         idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                         idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end-1)+ResNormt(4:3:end)+ResNormt(5:3:end)];
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end-1)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                                idx2(end)=[];
                             case 1
                                ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+[ResNormt(5:3:end) 0]];
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+[Rest(:,5:3:end) zeros(n,1)]];
                             case 2
                                ResNormt=[ResNormt(1)+ResNormt(2) ResNormt(3:3:end)+ResNormt(4:3:end)+ResNormt(5:3:end)];
                                Rest=[Rest(:,1)+Rest(:,2) Rest(:,3:3:end)+Rest(:,4:3:end)+Rest(:,5:3:end)];
                         end
                         ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                         Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);


                  elseif (k1==3 && k2==2) 
                         idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                         idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+ResNormt(3:3:end);
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+Rest(:,3:3:end);
                             case 1
                                ResNormt=ResNormt(1:3:end-1)+ResNormt(2:3:end)+ResNormt(3:3:end);
                                Rest=Rest(:,1:3:end-1)+Rest(:,2:3:end)+Rest(:,3:3:end);
                                idx2(end)=[];
                             case 2
                                ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
                                Rest=Rest(:,1:3:end)+Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n,1)];
                         end
                         ResNormz=ResNormz+sparse(idx1,idx2,ResNormt,nz,nz);
                         Resz=Resz+sparse(repelem(idx1,n,1),repelem(idx2,n,1),Rest(:),nz,nz);


                  elseif (k1==3 && k2==1) 
                         idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                         idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
                                Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n,1)];
                                idx2(1)=[];
                             case 1
                                ResNormt=ResNormt(2:3:end-2)+ResNormt(3:3:end)+ResNormt(4:3:end);
                                Rest=Rest(:,2:3:end-2)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                idx2(1)=[];
                             case 2
                                ResNormt=ResNormt(2:3:end-1)+ResNormt(3:3:end)+ResNormt(4:3:end);
                                Rest=Rest(:,2:3:end-1)+Rest(:,3:3:end)+Rest(:,4:3:end);
                                idx2(1)=[];
                                idx1(end)=[];
                         end
                         ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                         Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);


                  elseif (k1==6 && k2==1) 
                         idx1=idx(logical(data.idx_perturb_hes(:,k1)),i);
                         idx2=idx(logical(data.idx_perturb_hes(:,k2)),j);
                         switch mod(data.nps,3)
                             case 0
                                ResNormt=ResNormt(3:3:end)+[ResNormt(4:3:end) 0];
                                Rest=Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n,1)];
                                idx2(1)=[];
                             case {1,2}
                                ResNormt=ResNormt(3:3:end)+ResNormt(4:3:end);
                                Rest=Rest(:,3:3:end)+Rest(:,4:3:end);
                                idx2(1)=[];
                         end
                         ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                         Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);


                  end
              elseif (nt && j<=(np+nt) && k2==1) %k1=1 no pertub of X and U when i<=(np+nt)
                  if i<=(np+nt) 
                      if k1==1
                          ResNormt=sum(ResNormt,2);
                          Rest=sum(sum(Rest,2));
                          idx1=idx(1,j);
                          idx2=idx(1,i);
                          ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                          Resz=Resz+sparse(idx2,idx1,Rest,nz,nz);
                      end
                  else
                      if k1==1
                             if mod(data.nps,3)
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+ResNormt(4:3:end)];
                                Rest=[Rest(:,1) Rest(:,3:3:end)+Rest(:,4:3:end)];
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             else
                                ResNormt=[ResNormt(1) ResNormt(3:3:end)+[ResNormt(4:3:end) 0]];
                                Rest=[Rest(:,1) Rest(:,3:3:end)+[Rest(:,4:3:end) zeros(n,1)]];
                                idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                             end
                             idx1=idx(1,j)*ones(length(idx2),1);
                             ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                             Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      elseif k1==4
                            ResNormt=ResNormt(1:3:end);
                            Rest=Rest(:,1:3:end);
                            idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx1=idx(1,j)*ones(length(idx2),1);
                            ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                            Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      elseif k1==2
                            idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx1=idx(1,j)*ones(length(idx2),1);
                            switch mod(data.nps,3)
                                 case {0,2}
                                    ResNormt=ResNormt(1:3:end)+ResNormt(2:3:end);
                                    Rest=Rest(:,1:3:end)+Rest(:,2:3:end);
                                 case 1
                                    ResNormt=ResNormt(1:3:end)+[ResNormt(2:3:end) 0];
                                    Rest=Rest(:,1:3:end)+[Rest(:,2:3:end) zeros(n,1)];
                            end
                            ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                            Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      elseif k1==5
                            ResNormt=ResNormt(2:3:end);
                            Rest=Rest(:,2:3:end);
                            idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx1=idx(1,j)*ones(length(idx2),1);
                            ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                            Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      elseif k1==3
                            idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx1=idx(1,j)*ones(length(idx2),1);
                            switch mod(data.nps,3)
                                 case {0,1}
                                    ResNormt=ResNormt(2:3:end)+ResNormt(3:3:end);
                                    Rest=Rest(:,2:3:end)+Rest(:,3:3:end);
                                 case 2
                                    ResNormt=ResNormt(2:3:end)+[ResNormt(3:3:end) 0];
                                    Rest=Rest(:,2:3:end)+[Rest(:,3:3:end) zeros(n,1)];
                            end
                            ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                            Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      elseif k1==6
                            ResNormt=ResNormt(3:3:end);
                            Rest=Rest(:,3:3:end);
                            idx2=idx(logical(data.idx_perturb_hes(:,k1)),i);
                            idx1=idx(1,j)*ones(length(idx2),1);
                            ResNormz=ResNormz+sparse(idx2,idx1,ResNormt,nz,nz);
                            Resz=Resz+sparse(repelem(idx2,n,1),repelem(idx1,n,1),Rest(:),nz,nz);

                      end
                  end
              end
           end
        end
    end
end

% Return the Hessian of the Lagrangian
% -------------------------------------
hessc=gzz+bzz+Resz;
hessian=data.sigma*(Lzz+Ezz+ResNormz)+hessc;



%------------- END OF CODE --------------

