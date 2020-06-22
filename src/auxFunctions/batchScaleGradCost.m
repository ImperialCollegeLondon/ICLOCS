function [ dL,dE ] = batchScaleGradCost(dL,dE,data)
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
    if ~isempty(dL.dx)
        dL.dx=scale_variables_back( dL.dx, vdat.Xscale, 0 );
    end
    if ~isempty(dE.dx0)
        dE.dx0=scale_variables_back( dE.dx0, vdat.Xscale, 0 );
    end
    if ~isempty(dE.dxf)
        dE.dxf=scale_variables_back( dE.dxf, vdat.Xscale, 0 );
    end
    if ~isempty(dL.du)
        dL.du=scale_variables_back( dL.du, vdat.Uscale, 0);
    end
    if ~isempty(dE.du0)
        dE.du0=scale_variables_back( dE.du0, vdat.Uscale, 0);
    end
    if ~isempty(dE.duf)
        dE.duf=scale_variables_back( dE.duf, vdat.Uscale, 0);
    end

    if isfield(vdat,'Pscale')
        if ~isempty(dL.dp)
            dL.dp=scale_variables_back( dL.dp, vdat.Pscale, 0 );
        end
        if ~isempty(dE.dp)
            dE.dp=scale_variables_back( dE.dp, vdat.Pscale, 0 );
        end 
        
    end
end

end

