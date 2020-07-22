function [XU0f] = prepareForLinkageConstraintInfo(x0,xf,u0,uf,vdat)
%prepareForLinkageConstraintInfo - prepare the struct variable for linkage constraint
%
% Syntax:  [XU0f] = prepareForLinkageConstraintInfo(x0,xf,u0,uf,vdat)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

if isfield(vdat,'Xscale')
        XU0f.x0=scale_variables_back( x0', vdat.Xscale_back, vdat.Xshift )';
        XU0f.xf=scale_variables_back( xf', vdat.Xscale_back, vdat.Xshift )';
        XU0f.u0=scale_variables_back( u0', vdat.Uscale_back, vdat.Ushift )';
        XU0f.uf=scale_variables_back( uf', vdat.Uscale_back, vdat.Ushift )';
else
        XU0f.x0=x0;
        XU0f.xf=xf;
        XU0f.u0=u0;
        XU0f.uf=uf;
end




end

