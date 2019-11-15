function [ Ezz ] = hessian_LGR_AN_E( Ezz, HE, nz, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  idx=data.FD.index.Ey; 
  nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
  
    for i=1:nfd
        for j=1:i
          Ezz=Ezz+sparse(idx(i),idx(j),HE{j,i},nz,nz);
        end
    end



end

