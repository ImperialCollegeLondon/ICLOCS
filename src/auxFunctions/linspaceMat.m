function [ AB ] = linspaceMat( a,b,N )
%linspaceMat - Generates a matrix of linearly equally spaced points between column vectors a and b, and return in vector format
%
% Syntax:  [ AB ] = linspaceMat( a,b,N )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

if ~isempty(a) && ~isempty(b)
    dx = (b-a)/(N-1); 
    AB = repmat(dx,1,N); 
    AB(:,1) = a; 
    AB = cumsum(AB,2)';
    AB=AB(:);
else
    AB=[];
end

end

