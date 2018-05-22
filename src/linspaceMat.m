function [ AB ] = linspaceMat( A,B,N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(A) && ~isempty(B)
    dx = (B-A)/(N-1); 
    AB = repmat(dx,1,N); 
    AB(:,1) = A; 
    AB = cumsum(AB,2)';
    AB=AB(:);
else
    AB=[];
end

end

