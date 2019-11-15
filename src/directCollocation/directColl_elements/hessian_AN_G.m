function [ gzz ] = hessian_AN_G( gzz, Hg, M, nt, ng, nz, T, adjoint_g, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  idx=data.FD.index.g; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  
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


end

