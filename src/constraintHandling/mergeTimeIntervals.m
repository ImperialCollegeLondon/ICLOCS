function [ y ] = mergeTimeIntervals( x )
%mergeTimeIntervals - Merge overlapping time intervals
%
% Syntax:   [ y ] = mergeTimeIntervals( x )
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


x = sort(x,2);
nrow = size(x,1); 
[x,ind] = sort(x(:));
n = [(1:nrow) (1:nrow)]';
n = n(ind);
c = [ones(1,nrow) -ones(1,nrow)]';
c = c(ind);
csc = cumsum(c);
irit = find(csc==0);
ilef = [1; irit+1];
ilef(end) = []; 
y = [x(ilef) x(irit)];

end

