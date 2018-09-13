function [ U ] = computeU( z )
%computeU - obtain control trajectory from optimization variables
%
% Syntax:   [ U ] = computeU( z )
%
% Inputs:
%    z - optimization variables from NLP solution
%
% Outputs:
%    U - resultant trajectory for input variables
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------

[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});
U=reshape(data.map.Vu*z,m,N)';
U=U(1,:);
end

