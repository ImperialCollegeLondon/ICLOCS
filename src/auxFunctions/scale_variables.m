function [ scaledvari ] = scale_variables( vari, vscales, vshifts )
%scale_variables - scale the variables
%
% Syntax:  [ scaledvari ] = scale_variables( vari, vscales, vshifts )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
% 
%------------- BEGIN CODE --------------

scaledvari= bsxfun(@plus, bsxfun(@times, vari, vscales),vshifts);

end

