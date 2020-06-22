function [HL,HE,Hf,Hg,Hb] = batchScaleLagHessian(HL,HE,Hf,Hg,Hb,data)
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
    if ~isempty(HL)
        HL = cellfun(@times, HL, num2cell(1./vdat.Allscale_fgL_Mat), 'UniformOutput', false);
    end
    if ~isempty(Hf)
        Hf = cellfun(@times, Hf, num2cell(1./vdat.Allscale_fgL_Mat), 'UniformOutput', false);
    end
    if ~isempty(Hg)
        Hg = cellfun(@times, Hg, num2cell(1./vdat.Allscale_fgL_Mat), 'UniformOutput', false);
    end
    if ~isempty(HE)
        HE = cellfun(@times, HE, num2cell(1./vdat.Allscale_bE_Mat), 'UniformOutput', false);
    end
    if ~isempty(Hb)
        Hb = cellfun(@times, Hb, num2cell(1./vdat.Allscale_bE_Mat), 'UniformOutput', false);
    end
    
end

end

