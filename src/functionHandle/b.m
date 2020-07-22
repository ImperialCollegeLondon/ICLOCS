function bc=b(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
% Warp function
%------------- BEGIN CODE --------------
b_unscaled=vdat.functions_unscaled{6};
bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
if isfield(vdat,'Xscale')
    if ~isempty(bc)
        x0=scale_variables_back( x0', vdat.Xscale_back, vdat.Xshift )';
        xf=scale_variables_back( xf', vdat.Xscale_back, vdat.Xshift )';
        u0=scale_variables_back( u0', vdat.Uscale_back, vdat.Ushift )';
        uf=scale_variables_back( uf', vdat.Uscale_back, vdat.Ushift )';
        if isfield(vdat,'Pscale')
            p=scale_variables_back( p', vdat.Pscale_back, vdat.Pshift )';
        end
        bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
    end
end


%------------- END OF CODE ---------------------