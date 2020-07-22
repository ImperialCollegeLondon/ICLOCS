function [ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScale(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data)
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
    X=scale_variables( X, vdat.Xscale, vdat.Xshift );
    if ~isempty(Xr)
        Xr=scale_variables( Xr, vdat.Xscale, vdat.Xshift );
    end
    U=scale_variables( U, vdat.Uscale, vdat.Ushift );
    if ~isempty(Ur)
        Ur=scale_variables( Ur, vdat.Uscale, vdat.Ushift );
    end
    if isfield(vdat,'Pscale')
        P=scale_variables( P, vdat.Pscale, vdat.Pshift );
        p=scale_variables( p', vdat.Pscale, vdat.Pshift );
    end
    x0=scale_variables( x0', vdat.Xscale, vdat.Xshift );
    xf=scale_variables( xf', vdat.Xscale, vdat.Xshift );
    u0=scale_variables( u0', vdat.Uscale, vdat.Ushift );
    uf=scale_variables( uf', vdat.Uscale, vdat.Ushift );
end

end

