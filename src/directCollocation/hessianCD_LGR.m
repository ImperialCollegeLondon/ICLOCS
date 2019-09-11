function hessian=hessianCD_LGR(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%  It evaluates the Hessian of the Lagrangian with finite diferences
%  considering central difference formula
%
% Syntax: hessian=hessianCD_LGR(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation_LGR.m
%
% Outputs:
%    hessian - hessian information
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
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


% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,npd,npdu,npduidx,nps,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:17});
nrc=nrcl+nrcu+nrce;
nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
Gamma_1toN=reshape(data.lambda(1:M*n),M,n);
adjoint_f=Gamma_1toN;
adjoint_g=zeros(M,ng);
adjoint_g(logical(data.gActiveIdx))=data.lambda(n*M+1:n*M+ngActive);
adjoint_b=data.lambda((M*n+ngActive+nrc+1):end);

vdat=data.data;
fg=vdat.functionfg;

k0=tf+t0;
DT=tf-t0;
DTLP=repmat(data.t_segment_end,1,n);
DT_ratio_diff=repmat(data.tau_segment_ratio_diff,1,n);


e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size

% Compute fzz and gzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));
gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng && size(data.FD.index.f,2)==size(data.FD.index.g,2)
    idxf=data.FD.index.f;idxg=data.FD.index.g;nfd=size(idxf,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;
    
    for i=1:nfd
       for j=1:i
             if j==i
                 [dyn_1,gp1]=fg(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                 [dyn_2,gp2]=fg(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                 [dyn_o,go]=fg(X,U,P,DT/2*T+k0/2,vdat);
                 
                 fp1=(DTLP-et0(i)/2*DT_ratio_diff+etf(i)/2*DT_ratio_diff).*dyn_1;
                 fp2=(DTLP+et0(i)/2*DT_ratio_diff-etf(i)/2*DT_ratio_diff).*dyn_2;
                 fo=DTLP.*dyn_o;
                 ft=(fp2-2*fo+fp1)/e2;
                 fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 gt=(gp2-2*go+gp1).*adjoint_g/e2;
             else
                 [dyn_pp,gpp]=fg(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);
                 [dyn_pm,gpm]=fg(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                 [dyn_mm,gmm]=fg(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                 [dyn_mp,gmp]=fg(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                 
                 fpp=(DTLP+(-et0(i)-et0(j)+etf(i)+etf(j))/2*DT_ratio_diff).*dyn_pp;   
                 fpm=(DTLP+(-et0(i)+et0(j)+etf(i)-etf(j))/2*DT_ratio_diff).*dyn_pm; 
                 fmm=(DTLP+(+et0(i)+et0(j)-etf(i)-etf(j))/2*DT_ratio_diff).*dyn_mm;
                 fmp=(DTLP+(+et0(i)-et0(j)-etf(i)+etf(j))/2*DT_ratio_diff).*dyn_mp; 
                 ft=(fpp-fpm+fmm-fmp)/e2/4;
                 fzz=fzz+sparse(idxf(:,i),idxf(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 gt=(gpp-gpm+gmm-gmp).*adjoint_g/e2/4;
             end
             g_vect=reshape(gt,M*ng,1);
             gzz=gzz+sparse(idxg(data.gAllidx,i),idxg(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
       end
    end
    
else
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
             if j==i
                 fp1=(DTLP-et0(i)/2*DT_ratio_diff+etf(i)/2*DT_ratio_diff).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                 fp2=(DTLP+et0(i)/2*DT_ratio_diff-etf(i)/2*DT_ratio_diff).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                 fo=DTLP.*f(X,U,P,DT/2*T+k0/2,vdat);
                 ft=(fp2-2*fo+fp1)/e2;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             else
                 fpp=(DTLP+(-et0(i)-et0(j)+etf(i)+etf(j))/2*DT_ratio_diff).*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);   
                 fpm=(DTLP+(-et0(i)+et0(j)+etf(i)-etf(j))/2*DT_ratio_diff).*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                 fmm=(DTLP+(+et0(i)+et0(j)-etf(i)-etf(j))/2*DT_ratio_diff).*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                 fmp=(DTLP+(+et0(i)-et0(j)-etf(i)+etf(j))/2*DT_ratio_diff).*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                 ft=(fpp-fpm+fmm-fmp)/e2/4;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             end
       end
    end

    if ng
        idx=data.FD.index.g;nfd=size(idx,2);
        et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
        ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
        for i=1:nfd
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
             gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
            end
        end
    end
end
  
% Compute (w'L)zz
% ----------------

idx=data.FD.index.Ly;
nfd=size(idx,2);                               

et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

for i=1:nfd
dt0_1=e*et0{i};dtf_1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
    for j=1:i
      if j==i
        Lo=data.t_segment_end.*L(X,Xr,U,Ur,P,DT/2*T+k0/2,vdat);
        Lp1=(data.t_segment_end+dtf_1/2-dt0_1/2).*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf_1-dt0_1)/2*T+(k0+dtf_1+dt0_1)/2,vdat);
        Lp2=(data.t_segment_end-dtf_1/2+dt0_1/2).*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf_1+dt0_1)/2*T+(k0-dtf_1-dt0_1)/2,vdat);
        Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
      else
        dt0_2=e*et0{j};dtf_2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
        Lpp=(data.t_segment_end+dtf_1/2-dt0_1/2+dtf_2/2-dt0_2/2).*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,vdat);
        Lpm=(data.t_segment_end+dtf_1/2-dt0_1/2-dtf_2/2+dt0_2/2).*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,vdat);
        Lmp=(data.t_segment_end-dtf_1/2+dt0_1/2+dtf_2/2-dt0_2/2).*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,vdat);
        Lmm=(data.t_segment_end-dtf_1/2+dt0_1/2-dtf_2/2+dt0_2/2).*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,vdat);
        Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
      end
       Lzz=Lzz+sparse(idx(:,i),idx(:,j),Lt,nz,nz);
    end
end


% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
idx=data.FD.index.Ey;nfd=size(idx,2);                               
et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

for i=1:nfd
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


idx=data.FD.index.b;nfd=size(idx,2);
et0=e*data.FD.vector.b.et0;etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;et=e.*data.FD.vector.b.et;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;


for i=1:nfd
   for j=1:i
    if j==i
     bo=b(x0,xf,u0,uf,p,t0,tf,vdat,data.options,data.t_segment);
     bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat,data.options,data.t_segment+et(i,:)');
     bp2=b(x0-ex0(:,j),xf-exf(:,j),u0-eu0(:,j),uf-euf(:,j),p-ep(:,j),t0-et0(:,j),tf-etf(:,j),vdat,data.options,data.t_segment-et(j,:)');
     bt=(bp2-2*bo+bp1).*adjoint_b/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    else
    bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat,data.options,data.t_segment+et(i,:)'+et(j,:)');
    bpm=b(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
          uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat,data.options,data.t_segment+et(i,:)'-et(j,:)');
    bmp=b(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
          uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat,data.options,data.t_segment-et(i,:)'+et(j,:)');
    bmm=b(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
          uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat,data.options,data.t_segment-et(i,:)'-et(j,:)');
    bt=(bpp-bpm+bmm-bmp).*adjoint_b/e2/4; 
   
    bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    end
   end
end

end


% Return the Hessian of the Lagrangian
% -------------------------------------
hessian=data.sigma*(Lzz+Ezz)-fzz+gzz+bzz;


%------------- END OF CODE --------------

