function [ bzz ] = hessian_LGR_AN_B( bzz, Hb, nz, adjoint_b, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

idx=data.FD.index.b; 
nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
for i=1:nfd
   for j=1:i 
      bzz=bzz+sparse(idx(:,i),idx(:,j),Hb{j,i}.*adjoint_b,nz,nz);
   end
end


end

