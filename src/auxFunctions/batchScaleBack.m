function [ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data)
%batchScaleBack - Scale the variables back into original dimension in batches
%
% Syntax:  [ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data)
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


vdat=data.data;
if isfield(vdat,'Xscale')
    X=scale_variables_back( X, vdat.Xscale, vdat.Xshift );
    if ~isempty(Xr)
        Xr=scale_variables_back( Xr, vdat.Xscale, vdat.Xshift );
    end
    U=scale_variables_back( U, vdat.Uscale, vdat.Ushift );
    if ~isempty(Ur)
        Ur=scale_variables_back( Ur, vdat.Uscale, vdat.Ushift );
    end
    if isfield(vdat,'Pscale')
        P=scale_variables_back( P, vdat.Pscale, vdat.Pshift );
        p=scale_variables_back( p', vdat.Pscale, vdat.Pshift );
    end
    x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift );
    xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift );
    u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift );
    uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift );
end

end

