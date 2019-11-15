function [ Lzz ] = hessian_LGR_AN_wL( dL, Lzz, HL, M, nt, nz, alpha_j, beta_j, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 idx=data.FD.index.Ly; 
 nfd=size(idx,2);  
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
    
    
    


end

