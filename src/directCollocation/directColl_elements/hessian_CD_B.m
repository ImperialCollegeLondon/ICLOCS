function [ bzz ] = hessian_CD_B( bzz, nz, b, x0, xf, u0, uf, p, t0, tf, e, e2, adjoint, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx=data.FD.index.b;nfd=size(idx,2);
et0=e*data.FD.vector.b.et0;etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;%ez=e*data.FD.vector.b.ez;

for i=1:nfd
   for j=1:i
    if j==i
     bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
     bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),t0+et0(:,i),tf+etf(:,i),vdat);
     bp2=b(x0-ex0(:,j),xf-exf(:,j),u0-eu0(:,j),uf-euf(:,j),p-ep(:,j),t0-et0(:,i),tf-etf(:,j),vdat);
     bt=(bp2-2*bo+bp1).*adjoint'/e2; 
     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    else
    bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),t0+et0(:,i)+et0(:,j),tf+etf(:,i)+etf(:,j),vdat);
    bpm=b(x0+ex0(:,i)-ex0(:,j),xf+exf(:,i)-exf(:,j),u0+eu0(:,i)-eu0(:,j),...
          uf+euf(:,i)-euf(:,j),p+ep(:,i)-ep(:,j),t0+et0(:,i)-et0(:,j),tf+etf(:,i)-etf(:,j),vdat);
    bmp=b(x0-ex0(:,i)+ex0(:,j),xf-exf(:,i)+exf(:,j),u0-eu0(:,i)+eu0(:,j),...
          uf-euf(:,i)+euf(:,j),p-ep(:,i)+ep(:,j),t0-et0(:,i)+et0(:,j),tf-etf(:,i)+etf(:,j),vdat);
    bmm=b(x0-ex0(:,i)-ex0(:,j),xf-exf(:,i)-exf(:,j),u0-eu0(:,i)-eu0(:,j),...
          uf-euf(:,i)-euf(:,j),p-ep(:,i)-ep(:,j),t0-et0(:,i)-et0(:,j),tf-etf(:,i)-etf(:,j),vdat);
    bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 
   
    bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
    end
   end
end

% for i=1:nfd
%    for j=1:i
%     if j==i
%      bo=b(x0,xf,u0,uf,p,t0,tf,vdat);
%      bp1=b(x0+ex0(:,i)*ez(i),xf+exf(:,i)*ez(i),u0+eu0(:,i)*ez(i),uf+euf(:,i)*ez(i),p+ep(:,i)*ez(i),t0+et0(:,i)*ez(i),tf+etf(:,i)*ez(i),vdat);
%      bp2=b(x0-ex0(:,j)*ez(i),xf-exf(:,j)*ez(i),u0-eu0(:,j)*ez(i),uf-euf(:,j)*ez(i),p-ep(:,j)*ez(i),t0-et0(:,i)*ez(i),tf-etf(:,j)*ez(i),vdat);
%      bt=(bp2-2*bo+bp1).*adjoint'/e2; 
%      bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
%     else
%     bpp=b(x0+ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)+exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
%           uf+euf(:,i)*ez(i)+euf(:,j)*ez(j),p+ep(:,i)*ez(i)+ep(:,j)*ez(j),t0+et0(:,i)*ez(i)+et0(:,j)*ez(j),tf+etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
%     bpm=b(x0+ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf+exf(:,i)*ez(i)-exf(:,j)*ez(j),u0+eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
%           uf+euf(:,i)*ez(i)-euf(:,j)*ez(j),p+ep(:,i)*ez(i)-ep(:,j)*ez(j),t0+et0(:,i)*ez(i)-et0(:,j)*ez(j),tf+etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
%     bmp=b(x0-ex0(:,i)*ez(i)+ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)+exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)+eu0(:,j)*ez(j),...
%           uf-euf(:,i)*ez(i)+euf(:,j)*ez(j),p-ep(:,i)*ez(i)+ep(:,j)*ez(j),t0-et0(:,i)*ez(i)+et0(:,j)*ez(j),tf-etf(:,i)*ez(i)+etf(:,j)*ez(j),vdat);
%     bmm=b(x0-ex0(:,i)*ez(i)-ex0(:,j)*ez(j),xf-exf(:,i)*ez(i)-exf(:,j)*ez(j),u0-eu0(:,i)*ez(i)-eu0(:,j)*ez(j),...
%           uf-euf(:,i)*ez(i)-euf(:,j)*ez(j),p-ep(:,i)*ez(i)-ep(:,j)*ez(j),t0-et0(:,i)*ez(i)-et0(:,j)*ez(j),tf-etf(:,i)*ez(i)-etf(:,j)*ez(j),vdat);
%     bt=(bpp-bpm+bmm-bmp).*adjoint'/e2/4; 
%    
%     bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz); 
%     end
%    end
% end


end

