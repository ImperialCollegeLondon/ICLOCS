
function [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
% hessianAN_LGR - Return the Hessian of the Lagrangian when the analytic option has been selected
%
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianAN_LGR(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    Lzz - hessian of wL wrt z
%    Ezz - hessian of E wrt z
%    fzz - hessian of f wrt z
%    gzz - hessian of g wrt z
%    bzz - hessian of b wrt z
%
% Other m-files required: hessianLagrangian.m
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

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
t=(tf-t0)/2*T+(tf+t0)/2;
alpha=(1-T)/2;
beta=(1+T)/2;
hessianLagrangian=data.analyticDeriv.hessianLagrangian;



[HL,HE,Hf,Hg,Hb]=hessianLagrangian(X,U,P,t,E,x0,xf,u0,uf,p,t0,tf,data);



% Define some useful variables
[nt,np,n,m,ng,nb,M]=deal(data.sizes{1:7});
nz=nt+np+(M+1)*n+M*m;                              % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
Gamma_1toN=reshape(data.lambda(1:M*n),M,n);
adjoint_f=Gamma_1toN;
adjoint_g=reshape(data.lambda((M*n+1):M*(n+ng)),M,ng);
adjoint_b=data.lambda((M*(n+ng)+1):end);

vdat=data.data;
k0=tf+t0;
DT=tf-t0;
DTLP=repmat(data.t_segment_end,1,n);
alpha_j=kron(alpha,ones(1,n));
beta_j=kron(beta,ones(1,n));


% Compute fzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));  % Allocate some memory


if ~isempty(Hf)
  idx=data.FD.index.f; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  
  if ~isempty(df{end})
      for i=1:nfd
       for j=1:i
          if j==(nfd-nt+1) && length(df)==nfd
             ft=(df{i}+DTLP.*Hf{j,i}.*alpha_j).*adjoint_f;
          elseif j==(nfd) && length(df)==nfd
             ft=(df{i}+DTLP.*Hf{j,i}.*beta_j).*adjoint_f;
          else 
            ft=DTLP.*Hf{j,i}.*adjoint_f;
          end
          fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft,M*n,1),nz,nz);
        end
      end
  else
      for i=1:nfd-nt
       for j=1:i
          ft=DTLP.*Hf{j,i}.*adjoint_f;
          fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft,M*n,1),nz,nz);
        end
      end
      
        % numerical derivative w.r.t time
        et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
        ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

        for i=nfd-nt+1:nfd
           for j=1:i
                 if j==i
                     fp1=(DTLP-et0(i)/2+etf(i)/2).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                     fp2=(DTLP+et0(i)/2-etf(i)/2).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                     fo=DTLP.*f(X,U,P,DT/2*T+k0/2,vdat);
                     ft=(fp2-2*fo+fp1)/e2;
                     fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 else
                     fpp=(DTLP-et0(i)/2-et0(j)/2+etf(i)/2+etf(j)/2).*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);   
                     fpm=(DTLP-et0(i)/2+et0(j)/2+etf(i)/2-etf(j)/2).*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                     fmm=(DTLP+et0(i)/2+et0(j)/2-etf(i)/2-etf(j)/2).*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                     fmp=(DTLP+et0(i)/2-et0(j)/2-etf(i)/2+etf(j)/2).*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                     ft=(fpp-fpm+fmm-fmp)/e2/4;
                     fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
                 end
           end
        end
  end
  

else
% If Hf is empty the hessian related to the dynamical system is computed numerically    
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    for i=1:nfd
       for j=1:i
             if j==i
                 fp1=(DTLP-et0(i)/2+etf(i)/2).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)-et0(i))/2.*T+(k0+etf(i)+et0(i))/2,vdat);
                 fp2=(DTLP+et0(i)/2-etf(i)/2).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(DT-etf(i)+et0(i))/2.*T+(k0-etf(i)-et0(i))/2,vdat);
                 fo=DTLP.*f(X,U,P,DT/2*T+k0/2,vdat);
                 ft=(fp2-2*fo+fp1)/e2;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             else
                 fpp=(DTLP-et0(i)/2-et0(j)/2+etf(i)/2+etf(j)/2).*f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT-et0(i)-et0(j)+etf(i)+etf(j))/2.*T+(k0+et0(i)+et0(j)+etf(i)+etf(j))/2,vdat);   
                 fpm=(DTLP-et0(i)/2+et0(j)/2+etf(i)/2-etf(j)/2).*f(X+(ex{i}-ex{j})*e,U+(eu{i}-eu{j})*e,P+(ep{i}-ep{j})*e,(DT-et0(i)+et0(j)+etf(i)-etf(j))/2.*T+(k0+et0(i)-et0(j)+etf(i)-etf(j))/2,vdat); 
                 fmm=(DTLP+et0(i)/2+et0(j)/2-etf(i)/2-etf(j)/2).*f(X-(ex{i}+ex{j})*e,U-(eu{i}+eu{j})*e,P-(ep{i}+ep{j})*e,(DT+et0(i)+et0(j)-etf(i)-etf(j))/2.*T+(k0-et0(i)-et0(j)-etf(i)-etf(j))/2,vdat);
                 fmp=(DTLP+et0(i)/2-et0(j)/2-etf(i)/2+etf(j)/2).*f(X-(ex{i}-ex{j})*e,U-(eu{i}-eu{j})*e,P-(ep{i}-ep{j})*e,(DT+et0(i)-et0(j)-etf(i)+etf(j))/2.*T+(k0-et0(i)+et0(j)-etf(i)+etf(j))/2,vdat); 
                 ft=(fpp-fpm+fmm-fmp)/e2/4;
                 fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft.*adjoint_f,M*n,1),nz,nz);
             end
       end
    end
end



% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

% If there are path constraints (ng=1), their hessian is computed
% numerically whenever Hg is empty, otherwise the evaluation of the analytic espression 
% is exploited

if ng
  idx=data.FD.index.g; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  if ~isempty(Hg)
    [~,dg,~]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
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
  else    
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
         gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt,M*ng,1),nz,nz);
       end
    end
  end
end

% Compute (w'L)zz
% ----------------

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
idx=data.FD.index.Ly; 
nfd=size(idx,2);  
if (~isempty(HL))
  for i=1:(nfd-nt)
   for j=1:i
        if j==(nfd-nt+1)
            Lt=(dL{i}+data.t_segment_end.*HL{j,i}.*alpha_j).*data.map.w;
        elseif j==(nfd)   
            Lt=(dL{i}+data.t_segment_end.*HL{j,i}.*beta_j).*data.map.w;
        else    
            Lt=data.t_segment_end/2.*HL{j,i}.*data.map.w;
        end
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt,M,1),nz,nz);
    end
  end
    
else
    
% If HL is empty the hessian of the cost is computed numerically
    et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
    ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;
    for i=1:nfd
    dt0_1=e*et0{i};dtf_1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};
        for j=1:i
          if j==i;
            Lo=DTLP*L(X,Xr,U,Ur,P,DT/2*T+k0/2,vdat);
            Lp1=(DTLP+dtf_1-dt0_1)/2*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dtf_1-dt0_1)/2*T+(k0+dtf_1+dt0_1),vdat);
            Lp2=(DTLP-dtf_1+dt0_1)/2*L(X-dx1,Xr,U-du1,Ur,P-dp1,(DT-dtf_1+dt0_1)/2*T+(k0-dtf_1-dt0_1),vdat);
            Lt=(Lp2-2*Lo+Lp1).*data.map.w/e2;
          else
            dt0_2=e*et0{j};dtf_2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};
            Lpp=(DTLP+dtf_1-dt0_1+dtf_2-dt0_2)/2*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dtf_1-dt0_1+dtf_2-dt0_2)/2*T+(k0+dtf_1+dt0_1+dtf_2+dt0_2)/2,vdat);
            Lpm=(DTLP+dtf_1-dt0_1-dtf_2+dt0_2)/2*L(X+dx1-dx2,Xr,U+du1-du2,Ur,P+dp1-dp2,(DT+dtf_1-dt0_1-dtf_2+dt0_2)/2*T+(k0+dtf_1+dt0_1-dtf_2-dt0_2)/2,vdat);
            Lmp=(DTLP-dtf_1+dt0_1+dtf_2-dt0_2)/2*L(X-dx1+dx2,Xr,U-du1+du2,Ur,P-dp1+dp2,(DT-dtf_1+dt0_1+dtf_2-dt0_2)/2*T+(k0-dtf_1-dt0_1+dtf_2+dt0_2)/2,vdat);
            Lmm=(DTLP-dtf_1+dt0_1-dtf_2+dt0_2)/2*L(X-dx1-dx2,Xr,U-du1-du2,Ur,P-dp1-dp2,(DT-dtf_1+dt0_1-dtf_2+dt0_2)/2*T+(k0-dtf_1-dt0_1-dtf_2-dt0_2)/2,vdat);
            Lt=(Lpp-Lmp+Lmm-Lpm).*data.map.w/e2/4;
          end
           Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt,M,1),nz,nz);
        end
    end
end



% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));

if ~isempty(HE)
  idx=data.FD.index.Ey; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  
    for i=1:nfd
        for j=1:i
          Ezz=Ezz+sparse(idx(i),idx(j),HE{j,i},nz,nz);
        end
    end

else    
idx=data.FD.index.Ey;nfd=size(idx,2);                               
et0=e*data.FD.vector.Ey.et0;etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;
Eo=E(x0,xf,u0,uf,p,t0,tf,vdat);
for i=1:nfd
    Ep1=E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat);

    for j=1:i
    
    if j==i;Ep2=Ep1;else
    Ep2=E(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),t0+et0(:,j),tf+etf(:,j),vdat);
    end

    Epp=E(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
    Ezz=Ezz+sparse(idx(i),idx(j),(Epp+Eo-Ep1-Ep2)/e2,nz,nz);
    end
end

end


% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));


if nb
  if ~isempty(Hb)
       idx=data.FD.index.b; 
       nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
       for i=1:nfd
           for j=1:i 
              bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint_b,nz,nz);
           end
       end
  else   
      idx=data.FD.index.b;nfd=size(idx,2);
      et0=e*data.FD.vector.b.et0;etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
      ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
      exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;
      bo=b(x0,xf,u0,uf,p,tf,vdat);
      for i=1:nfd
        bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat);
        for j=1:i
         if j==i;bp2=bp1;else
         bp2=b(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),t0+et0(:,j),tf+etf(:,j),vdat);
         end
         bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
              uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
         bt=(bpp-bp2+bo-bp1).*adjoint_b/e2; 
         bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz);
        end
      end

 end
end





