function [ x,u,p ] = scale_Problem_back( x,u,p )
%scale_Problem_back - scale the problem back
%
% Syntax:  [ x,u,p ] = scale_Problem_back( x,u,p )
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

if isfield(data.dataNLP.data,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(data.dataNLP.data,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
end
        
end

