
function hessian=resMin_hessianCD_LGR_ModeMinCost_Adigator(L,g,X_Np1,U_Np1,P,T_Np1,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment,data,Res_vec)
%resMin_hessianCD_LGR_ModeMinCost_Adigator - Hessian computation for
%integrated residual minimization (alternating method: cost
%minimization) with Adigator
%
% Syntax:   hessian=resMin_hessianCD_LGR_ModeMinCost_Adigator(L,g,X_Np1,U_Np1,P,T_Np1,E,b,x0,xf,u0,uf,p,t0,tf,t_segment_end,t_segment,data,Res_vec)
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
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(dataNLP.sizes{1:17});
nrc=nrcl+nrcu+nrce;
if data.free_time
    nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
else
    nz=np+(M+1)*n+M*m;                              % Length of the primal variable
end
% 

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


    if free_time
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
adjoint_Res=lambda(end-n+1:end);
adjoint_Res_vec=adjoint_Res(Res_vec.dYdY_location(:,1));
Resz=sparse(Res_vec.dYdY_location(:,2),Res_vec.dYdY_location(:,3),Res_vec.dYdY.*adjoint_Res_vec,Res_vec.dYdY_size(2),Res_vec.dYdY_size(3));


% Return the Hessian of the Lagrangian
% -------------------------------------?
hessc=gzz+bzz+Resz;
hessian=sigma*(Lzz+Ezz)+hessc;

