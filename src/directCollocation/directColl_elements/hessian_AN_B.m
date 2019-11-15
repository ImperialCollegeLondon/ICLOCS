function [ bzz ] = hessian_AN_B( bzz, Hb, nz, nt, adjoint, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%   adjoint=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';  

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


end

