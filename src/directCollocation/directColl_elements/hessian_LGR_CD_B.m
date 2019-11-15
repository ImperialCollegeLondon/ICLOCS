function [ bzz ] = hessian_LGR_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint_b, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

