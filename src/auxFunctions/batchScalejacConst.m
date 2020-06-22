function [ df,dg,db ] = batchScalejacConst(df,dg,db,data)
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
    if df.flag 
        if ~isempty(df.dx)
            for i=1:length(df.dx)
                df.dx{i}=df.dx{i}.*vdat.Xscale;
            end
            df.dx = cellfun(@times, df.dx, num2cell(1./vdat.Xscale), 'UniformOutput', false);
        end
        if ~isempty(df.du)
            for i=1:length(df.du)
                df.du{i}=df.du{i}.*vdat.Xscale;
            end
            df.du = cellfun(@times, df.du, num2cell(1./vdat.Uscale), 'UniformOutput', false);
        end
        if isfield(vdat,'Pscale')
            if ~isempty(df.dp)
                df.dp = cellfun(@times, df.dp, num2cell(1./vdat.Pscale), 'UniformOutput', false);
            end
        end
    end
    if dg.flag 
        if ~isempty(dg.dx)
            dg.dx = cellfun(@times, dg.dx, num2cell(1./vdat.Xscale), 'UniformOutput', false);
        end
        if ~isempty(dg.du)
            dg.du = cellfun(@times, dg.du, num2cell(1./vdat.Uscale), 'UniformOutput', false);
        end
        if isfield(vdat,'Pscale')
            if ~isempty(dg.dp)
                dg.dp = cellfun(@times, dg.dp, num2cell(1./vdat.Pscale), 'UniformOutput', false);
            end
        end
    end
    if db.flag 
        if ~isempty(db.dx0)
            db.dx0 = db.dx0./vdat.Xscale;
        end
        if ~isempty(db.dxf)
            db.dxf = db.dxf./vdat.Xscale;
        end
        if ~isempty(db.du0)
            db.du0 = db.du0./vdat.Uscale;
        end
        if ~isempty(db.duf)
            db.duf = db.duf./vdat.Uscale;
        end
        
        if isfield(vdat,'Pscale')
            if ~isempty(db.dp)
                db.dp = db.dp./vdat.Pscale;
            end
        end
    end
end

end

