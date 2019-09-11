
function [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,df,dL,X,U,P,T,E,b,x0,xf,u0,uf,p,t0,tf,data)
                               
% hessianAN - Return the Hessian of the Lagrangian when the
%             analytic option has been selected
%
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)
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
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk


%------------- BEGIN CODE ---------------


% Define some useful variables

e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size
t=(tf-t0)*T+t0;
hessianLagrangian=data.analyticDeriv.hessianLagrangian;



[HL,HE,Hf,Hg,Hb]=hessianLagrangian(X,U,P,t,E,x0,xf,u0,uf,p,t0,tf,data);



% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns,nrcl,nrcu,nrce,ngActive]=deal(data.sizes{1:13});
nz=nt+np+M*n+N*m;                           % Length of the primal variable


Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';

vdat=data.data;
DT=tf-t0;
Tj=kron(T,ones(1,n));
% Compute fzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));  % Allocate some memory


if ~isempty(Hf)
      idx=data.FD.index.f; 
      nfd=size(idx,2);
      
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
else
% If Hf is empty the hessian related to the dynamical system is computed numerically    
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

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



% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

% If there are path constraints (ng=1), their hessian is computed
% numerically whenever Hg is empty, otherwise the evaluation of the analytic espression 
% is exploited

if ng
    adjoint_g=zeros(ng,M);
    adjoint_g(logical(data.gActiveIdx'))=lambda(n*M+1:n*M+ngActive);
    adjoint_g=adjoint_g';
  idx=data.FD.index.g; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  if ~isempty(Hg)
  Tj=kron(T,ones(1,ng));
    if nt==1
      gt=(Hg{1,1}).*(Tj.^2).*adjoint_g;  
      g_vect=reshape(gt',M*ng,1);
      gzz=gzz+sparse(idx(data.gAllidx,1),idx(data.gAllidx,1),g_vect(data.gAllidx),nz,nz);
    elseif nt==2
      gt0t0=-(Hg{1,1}).*(Tj.^2).*adjoint_g;  
      gtftf=(Hg{2,2}).*(Tj.^2).*adjoint_g;  
      gt0tf=-(Hg{1,2}).*(Tj.^2).*adjoint_g;  
      gt0t0_vect=reshape(gt0t0',M*ng,1);
      gtftf_vect=reshape(gtftf',M*ng,1);
      gt0tf_vect=reshape(gt0tf',M*ng,1);
      gzz=gzz+sparse(idx(data.gAllidx,1),idx(data.gAllidx,1),gt0t0_vect(data.gAllidx),nz,nz)+sparse(idx(data.gAllidx,2),idx(data.gAllidx,2),gtftf_vect(data.gAllidx),nz,nz)+sparse(idx(data.gAllidx,2),idx(data.gAllidx,1),gt0tf_vect(data.gAllidx),nz,nz);
    end    
    for i=1+nt:nfd
    for j=1:i
     if (nt==1&&(j==1))
      gt=(Hg{j,i}).*(Tj).*adjoint_g;  
     elseif (nt==2&&(j<=nt))
        if j==1
            gt=-(Hg{j,i}).*(Tj).*adjoint_g;  
        elseif j==2
            gt=(Hg{j,i}).*(Tj).*adjoint_g;  
        end
     else   
      gt=Hg{j,i}.*adjoint_g;
     end
     g_vect=reshape(gt',M*ng,1);
     gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
    end
   end
  else    
    idx=data.FD.index.g;nfd=size(idx,2);
    et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
    ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;ez=data.FD.vector.g.ez;
    for i=1:nfd
       for j=1:i
        if j==i
         go=g(X,U,P,DT*T+t0,vdat);
         gp1=g(X+ex{i}*e*ez(i),U+eu{i}*e*ez(i),P+ep{i}*e*ez(i),(DT+etf(i)-et0(i)).*T+t0+et0(i),vdat);
         gp2=g(X-ex{j}*e*ez(i),U-eu{j}*e*ez(i),P-ep{j}*e*ez(i),(DT-etf(j)+et0(j)).*T+t0-et0(j),vdat);
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
         gzz=gzz+sparse(idx(data.gAllidx,i),idx(data.gAllidx,j),g_vect(data.gAllidx),nz,nz);
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
 if nfd&&nt==1
   Lt=(2*dL{1}.*T+DT*HL{1,1}.*(T.^2)).*data.map.w;  
   Lzz=Lzz+sparse(idx(:,1),idx(:,1),reshape(Lt',M,1),nz,nz);  
 elseif nfd&&nt==2
   Lt0t0=-(2*dL{1}.*T+DT*HL{1,1}.*(T.^2)).*data.map.w;  
   Ltftf=(2*dL{2}.*T+DT*HL{2,2}.*(T.^2)).*data.map.w;  
   Lt0tf=-(2*dL{1}.*T+DT*HL{1,2}.*(T.^2)).*data.map.w;  
   Lzz=Lzz+sparse(idx(:,1),idx(:,1),reshape(Lt0t0',M,1),nz,nz)+sparse(idx(:,2),idx(:,2),reshape(Ltftf',M,1),nz,nz)+sparse(idx(:,2),idx(:,1),reshape(Lt0tf',M,1),nz,nz);  
  end    
  for i=1+nt:nfd
   for j=1:i
      if (nt==1&&(j==1))
        Lt=(dL{i}+DT*HL{j,i}.*T).*data.map.w;
      elseif (nt==2&&(j<=nt))
        if j==1
            Lt=-(dL{i}+DT*HL{j,i}.*T).*data.map.w;
        elseif j==2
            Lt=(dL{i}+DT*HL{j,i}.*T).*data.map.w;
        end
      else   
        Lt=DT*HL{j,i}.*data.map.w;
      end
     Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
    end
  end
else
% If HL is empty the hessian of the cost is computed numerically
                         
    et0=data.FD.vector.Ly.et0;etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
    ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;ez=data.FD.vector.Ly.ez;

    Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);

    for i=1:nfd
    dt01=e*et0{i};dtf1=e*etf{i};dp1=e*ep{i}*ez(i);dx1=e*ex{i}*ez(i); du1=e*eu{i}*ez(i);
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





% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));




if ~isempty(HE)
  idx=data.FD.index.Ey; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  if nfd&& (nt==1 || length(idx)==1)
   Ezz(idx(1),idx(1))=HE{1,1};
 elseif nfd&&nt==2
   Ezz(idx(1),idx(1))=-HE{1,1};
   Ezz(idx(2),idx(2))=HE{2,2};
   Ezz(idx(2),idx(1))=-HE{1,2};
  end    
for i=1+nt:nfd
    for j=1:i
        if (nt==1&&(j==1))
            Ezz(idx(i),idx(j))=HE{j,i};
        elseif (nt==2&&(j<=nt))
            if j==1
                Ezz(idx(i),idx(j))=-HE{j,i};
            elseif j==2
                Ezz(idx(i),idx(j))=HE{j,i};
            end
        else   
            Ezz(idx(i),idx(j))=HE{j,i};
        end
    end
end

else    
    idx=data.FD.index.Ey;nfd=size(idx,2);   
    et0=data.FD.vector.Ey.et0;etf=data.FD.vector.Ey.etf;ep=data.FD.vector.Ey.ep;
    ex0=data.FD.vector.Ey.ex0;eu0=data.FD.vector.Ey.eu0;
    exf=data.FD.vector.Ey.exf;euf=data.FD.vector.Ey.euf;ez=e*data.FD.vector.Ey.ez;

    for i=1:nfd
       for j=1:i
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

end


% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));


if nb
  adjoint=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';  
  if ~isempty(Hb)
   idx=data.FD.index.b; 
   nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
   
  if nfd&&nt==1
   bzz=bzz+sparse(idx(:,1),idx(:,1),Hb{1,1}.*adjoint',nz,nz);
 elseif nfd&&nt==2
   bzz=bzz+sparse(idx(:,1),idx(:,1),-Hb{1,1}.*adjoint',nz,nz)+sparse(idx(:,2),idx(:,2),Hb{2,2}.*adjoint',nz,nz)+sparse(idx(:,2),idx(:,1),-Hb{1,2}.*adjoint',nz,nz);
  end    
   for i=1+nt:nfd
    for j=1:i 
        if (nt==1&&(j==1))
            bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint',nz,nz);
        elseif (nt==2&&(j<=nt))
            if j==1
                bzz=bzz+sparse(idx(:,i),idx(:,j),-Hb{j,i}.*adjoint',nz,nz);
            elseif j==2
                bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint',nz,nz);
            end
        else   
            bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint',nz,nz);
        end
        
    end
   end
  else   
    idx=data.FD.index.b;nfd=size(idx,2);
    et0=data.FD.vector.b.et0;etf=data.FD.vector.b.etf;ep=data.FD.vector.b.ep;
    ex0=data.FD.vector.b.ex0;eu0=data.FD.vector.b.eu0;
    exf=data.FD.vector.b.exf;euf=data.FD.vector.b.euf;ez=e*data.FD.vector.b.ez;

    adjoint=data.lambda(n*M+ngActive+nrc+(~~nb):n*M+ngActive+nrc+nb).';



    for i=1:nfd
       for j=1:i
        if j==i
         bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
         bp1=b(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i)*ez(i),tf+etf(:,i)*ez(i),vdat);
         bp2=b(x0-ex0(:,j)*ez(i),xf-exf(:,j)*ez(i),u0-eu0(:,j)*ez(i),uf-euf(:,j)*ez(i),p-ep(:,j)*ez(i),t0-et0(:,i)*ez(i),tf-etf(:,j)*ez(i),vdat);
         bt=(bp2-2*bo+bp1).*adjoint'/e2; 
         bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        else
        bpp=b(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)+exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
              uf+euf(:,i)*ez(i)+euf(:,j)*ez(j),p+ep(:,i)*ez(i)+ep(:,j)*ez(j),t0+et0(:,i)*ez(i)+et0(:,j)*ez(j),tf+etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
        bpm=b(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)-exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
              uf+euf(:,i)*ez(i)-euf(:,j)*ez(j),p+ep(:,i)*ez(i)-ep(:,j)*ez(j),t0+et0(:,i)*ez(i)-et0(:,j)*ez(j),tf+etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
        bmp=b(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)+exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
              uf-euf(:,i)*ez(i)+euf(:,j)*ez(j),p-ep(:,i)*ez(i)+ep(:,j)*ez(j),t0-et0(:,i)*ez(i)+et0(:,j)*ez(j),tf-etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
        bmm=b(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)-exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
              uf-euf(:,i)*ez(i)-euf(:,j)*ez(j),p-ep(:,i)*ez(i)-ep(:,j)*ez(j),t0-et0(:,i)*ez(i)-et0(:,j)*ez(j),tf-etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
        bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 

        bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
        end
       end
    end


 end
end





